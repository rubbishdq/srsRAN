/*
 * Copyright 2013-2020 Software Radio Systems Limited
 *
 * This file is part of srsLTE.
 *
 * srsLTE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * srsLTE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * A copy of the GNU Affero General Public License can be found in
 * the LICENSE file in the top-level directory of this distribution
 * and at http://www.gnu.org/licenses/.
 *
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <fcntl.h>

#include "srslte/config.h"
#include "srslte/phy/ch_estimation/chest_ul.h"
#include "srslte/phy/dft/dft_precoding.h"
#include "srslte/phy/utils/convolution.h"
#include "srslte/phy/utils/vector.h"
#include "srslte/srslte.h"

#define NOF_REFS_SYM (q->cell.nof_prb * SRSLTE_NRE)
#define NOF_REFS_SF (NOF_REFS_SYM * 2) // 2 reference symbols per subframe

#define MAX_REFS_SYM (max_prb * SRSLTE_NRE)
#define MAX_REFS_SF (max_prb * SRSLTE_NRE * 2) // 2 reference symbols per subframe

/** 3GPP LTE Downlink channel estimator and equalizer.
 * Estimates the channel in the resource elements transmitting references and interpolates for the rest
 * of the resource grid.
 *
 * The equalizer uses the channel estimates to produce an estimation of the transmitted symbol.
 *
 * This object depends on the srslte_refsignal_t object for creating the LTE CSR signal.
 */

int srslte_chest_ul_init(srslte_chest_ul_t* q, uint32_t max_prb)
{
  int ret = SRSLTE_ERROR_INVALID_INPUTS;
  if (q != NULL) {
    bzero(q, sizeof(srslte_chest_ul_t));

    ret = srslte_refsignal_ul_init(&q->dmrs_signal, max_prb);
    if (ret != SRSLTE_SUCCESS) {
      ERROR("Error initializing CSR signal (%d)\n", ret);
      goto clean_exit;
    }

    q->tmp_noise = srslte_vec_cf_malloc(MAX_REFS_SF);
    if (!q->tmp_noise) {
      perror("malloc");
      goto clean_exit;
    }
    q->pilot_estimates = srslte_vec_cf_malloc(MAX_REFS_SF);
    if (!q->pilot_estimates) {
      perror("malloc");
      goto clean_exit;
    }
    for (int i = 0; i < 4; i++) {
      q->pilot_estimates_tmp[i] = srslte_vec_cf_malloc(MAX_REFS_SF);
      if (!q->pilot_estimates_tmp[i]) {
        perror("malloc");
        goto clean_exit;
      }
    }
    q->pilot_recv_signal = srslte_vec_cf_malloc(MAX_REFS_SF + 1);
    if (!q->pilot_recv_signal) {
      perror("malloc");
      goto clean_exit;
    }

    q->pilot_known_signal = srslte_vec_cf_malloc(MAX_REFS_SF + 1);
    if (!q->pilot_known_signal) {
      perror("malloc");
      goto clean_exit;
    }

    if (srslte_interp_linear_vector_init(&q->srslte_interp_linvec, MAX_REFS_SYM)) {
      ERROR("Error initializing vector interpolator\n");
      goto clean_exit;
    }

    q->smooth_filter_len = 3;
    srslte_chest_set_smooth_filter3_coeff(q->smooth_filter, 0.3333);

    q->dmrs_signal_configured = false;

    if (srslte_refsignal_dmrs_pusch_pregen_init(&q->dmrs_pregen, max_prb)) {
      ERROR("Error allocating memory for pregenerated signals\n");
      goto clean_exit;
    }
  }

  ret = SRSLTE_SUCCESS;

clean_exit:
  if (ret != SRSLTE_SUCCESS) {
    srslte_chest_ul_free(q);
  }
  return ret;
}

void srslte_chest_ul_free(srslte_chest_ul_t* q)
{
  srslte_refsignal_dmrs_pusch_pregen_free(&q->dmrs_signal, &q->dmrs_pregen);

  srslte_refsignal_ul_free(&q->dmrs_signal);
  if (q->tmp_noise) {
    free(q->tmp_noise);
  }
  srslte_interp_linear_vector_free(&q->srslte_interp_linvec);

  if (q->pilot_estimates) {
    free(q->pilot_estimates);
  }
  for (int i = 0; i < 4; i++) {
    if (q->pilot_estimates_tmp[i]) {
      free(q->pilot_estimates_tmp[i]);
    }
  }
  if (q->pilot_recv_signal) {
    free(q->pilot_recv_signal);
  }
  if (q->pilot_known_signal) {
    free(q->pilot_known_signal);
  }
  bzero(q, sizeof(srslte_chest_ul_t));
}

int srslte_chest_ul_res_init(srslte_chest_ul_res_t* q, uint32_t max_prb)
{
  bzero(q, sizeof(srslte_chest_ul_res_t));
  q->nof_re = SRSLTE_SF_LEN_RE(max_prb, SRSLTE_CP_NORM);
  q->ce     = srslte_vec_cf_malloc(q->nof_re);
  if (!q->ce) {
    perror("malloc");
    return -1;
  }
  return 0;
}

void srslte_chest_ul_res_set_identity(srslte_chest_ul_res_t* q)
{
  for (uint32_t i = 0; i < q->nof_re; i++) {
    q->ce[i] = 1.0;
  }
}

void srslte_chest_ul_res_free(srslte_chest_ul_res_t* q)
{
  if (q->ce) {
    free(q->ce);
  }
}

int srslte_chest_ul_set_cell(srslte_chest_ul_t* q, srslte_cell_t cell)
{
  int ret = SRSLTE_ERROR_INVALID_INPUTS;
  if (q != NULL && srslte_cell_isvalid(&cell)) {
    if (cell.id != q->cell.id || q->cell.nof_prb == 0) {
      q->cell = cell;
      ret     = srslte_refsignal_ul_set_cell(&q->dmrs_signal, cell);
      if (ret != SRSLTE_SUCCESS) {
        ERROR("Error initializing CSR signal (%d)\n", ret);
        return SRSLTE_ERROR;
      }

      if (srslte_interp_linear_vector_resize(&q->srslte_interp_linvec, NOF_REFS_SYM)) {
        ERROR("Error initializing vector interpolator\n");
        return SRSLTE_ERROR;
      }
    }
    ret = SRSLTE_SUCCESS;
  }
  return ret;
}

void srslte_chest_ul_pregen(srslte_chest_ul_t*                 q,
                            srslte_refsignal_dmrs_pusch_cfg_t* cfg,
                            srslte_refsignal_srs_cfg_t*        srs_cfg)
{
  srslte_refsignal_dmrs_pusch_pregen(&q->dmrs_signal, &q->dmrs_pregen, cfg);
  q->dmrs_signal_configured = true;

  if (srs_cfg) {
    srslte_refsignal_srs_pregen(&q->dmrs_signal, &q->srs_pregen, srs_cfg, cfg);
    q->srs_signal_configured = true;
  }
}

/* Uses the difference between the averaged and non-averaged pilot estimates */
static float estimate_noise_pilots(srslte_chest_ul_t* q, cf_t* ce, uint32_t nslots, uint32_t nrefs, uint32_t n_prb[2])
{

  float power = 0;
  for (int i = 0; i < nslots; i++) {
    power += srslte_chest_estimate_noise_pilots(
        &q->pilot_estimates[i * nrefs],
        &ce[SRSLTE_REFSIGNAL_UL_L(i, q->cell.cp) * q->cell.nof_prb * SRSLTE_NRE + n_prb[i] * SRSLTE_NRE],
        q->tmp_noise,
        nrefs);
  }

  power /= nslots;

  if (q->smooth_filter_len == 3) {
    // Calibrated for filter length 3
    float w = q->smooth_filter[0];
    float a = 7.419 * w * w + 0.1117 * w - 0.005387;
    return (power / (a * 0.8));
  } else {
    return power;
  }
}

// The interpolator currently only supports same frequency allocation for each subframe
#define cesymb(i) ce[SRSLTE_RE_IDX(q->cell.nof_prb, i, n_prb[0] * SRSLTE_NRE)]
static void interpolate_pilots(srslte_chest_ul_t* q, cf_t* ce, uint32_t nslots, uint32_t nrefs, uint32_t n_prb[2])
{
#ifdef DO_LINEAR_INTERPOLATION
  uint32_t L1 = SRSLTE_REFSIGNAL_UL_L(0, q->cell.cp);
  uint32_t L2 = SRSLTE_REFSIGNAL_UL_L(1, q->cell.cp);
  uint32_t NL = 2 * SRSLTE_CP_NSYMB(q->cell.cp);

  /* Interpolate in the time domain between symbols */
  srslte_interp_linear_vector3(
      &q->srslte_interp_linvec, &cesymb(L2), &cesymb(L1), &cesymb(L1), &cesymb(L1 - 1), (L2 - L1), L1, false, nrefs);
  srslte_interp_linear_vector3(
      &q->srslte_interp_linvec, &cesymb(L1), &cesymb(L2), NULL, &cesymb(L1 + 1), (L2 - L1), (L2 - L1) - 1, true, nrefs);
  srslte_interp_linear_vector3(&q->srslte_interp_linvec,
                               &cesymb(L1),
                               &cesymb(L2),
                               &cesymb(L2),
                               &cesymb(L2 + 1),
                               (L2 - L1),
                               (NL - L2) - 1,
                               true,
                               nrefs);
#else
  // Instead of a linear interpolation, we just copy the estimates to all symbols in that subframe
  for (int s = 0; s < nslots; s++) {
    for (int i = 0; i < SRSLTE_CP_NSYMB(q->cell.cp); i++) {
      int src_symb = SRSLTE_REFSIGNAL_UL_L(s, q->cell.cp);
      int dst_symb = i + s * SRSLTE_CP_NSYMB(q->cell.cp);

      // skip the symbol with the estimates
      if (dst_symb != src_symb) {
        srslte_vec_cf_copy(&ce[(dst_symb * q->cell.nof_prb + n_prb[s]) * SRSLTE_NRE],
                           &ce[(src_symb * q->cell.nof_prb + n_prb[s]) * SRSLTE_NRE],
                           nrefs);
      }
    }
  }
#endif
}

static void
average_pilots(srslte_chest_ul_t* q, cf_t* input, cf_t* ce, uint32_t nslots, uint32_t nrefs, uint32_t n_prb[2])
{
  for (uint32_t i = 0; i < nslots; i++) {
    srslte_chest_average_pilots(
        &input[i * nrefs],
        &ce[SRSLTE_REFSIGNAL_UL_L(i, q->cell.cp) * q->cell.nof_prb * SRSLTE_NRE + n_prb[i] * SRSLTE_NRE],
        q->smooth_filter,
        nrefs,
        1,
        q->smooth_filter_len);
  }
}

/**
 * Generic PUSCH and DMRS channel estimation. It assumes q->pilot_estimates has been populated with the Least Square
 * Estimates
 *
 * @param q Uplink Channel estimation instance
 * @param nslots number of slots (2 for DMRS, 1 for SRS)
 * @param nrefs_sym number of reference resource elements per symbols (depends on configuration)
 * @param stride sub-carrier distance between reference signal resource elements (1 for DMRS, 2 for SRS)
 * @param meas_ta_en enables or disables the Time Alignment error measurement
 * @param write_estimates Write channel estimation in res, (true for DMRS and false for SRS)
 * @param n_prb Resource block start for the grant, set to zero for Sounding Reference Signals
 * @param res UL channel estimation result
 */
static void chest_ul_estimate(srslte_chest_ul_t*     q,
                              uint32_t               nslots,
                              uint32_t               nrefs_sym,
                              uint32_t               stride,
                              bool                   meas_ta_en,
                              bool                   write_estimates,
                              uint32_t               n_prb[SRSLTE_NOF_SLOTS_PER_SF],
                              srslte_chest_ul_res_t* res)
{
  // Calculate time alignment error
  float ta_err = 0.0f;
  if (meas_ta_en) {
    for (int i = 0; i < nslots; i++) {
      ta_err += srslte_vec_estimate_frequency(&q->pilot_estimates[i * nrefs_sym], nrefs_sym) / nslots;
    }
  }

  // Calculate actual time alignment error in micro-seconds
  if (isnormal(ta_err) && stride > 0) {
    ta_err /= (float)stride;                     // Divide by the pilot spacing
    ta_err /= 15e3f;                             // Convert from normalized frequency to seconds
    ta_err *= 1e6f;                              // Convert to micro-seconds
    ta_err     = roundf(ta_err * 10.0f) / 10.0f; // Round to one tenth of micro-second
    res->ta_us = ta_err;
  } else {
    res->ta_us = 0.0f;
  }

  //printf("%6f\n", res->ta_us);

  // Check if intra-subframe frequency hopping is enabled
  if (n_prb[0] != n_prb[1]) {
    ERROR("ERROR: intra-subframe frequency hopping not supported in the estimator!!\n");
  }

  if (res->ce != NULL) {
    if (q->smooth_filter_len > 0) {
      average_pilots(q, q->pilot_estimates, res->ce, nslots, nrefs_sym, n_prb);

      if (write_estimates) {
        interpolate_pilots(q, res->ce, nslots, nrefs_sym, n_prb);
      }

      // If averaging, compute noise from difference between received and averaged estimates
      res->noise_estimate = estimate_noise_pilots(q, res->ce, nslots, nrefs_sym, n_prb);
    } else {
      // Copy estimates to CE vector without averaging
      for (int i = 0; i < nslots; i++) {
        srslte_vec_cf_copy(
            &res->ce[SRSLTE_REFSIGNAL_UL_L(i, q->cell.cp) * q->cell.nof_prb * SRSLTE_NRE + n_prb[i] * SRSLTE_NRE],
            &q->pilot_estimates[i * nrefs_sym],
            nrefs_sym);
      }
      if (write_estimates) {
        interpolate_pilots(q, res->ce, nslots, nrefs_sym, n_prb);
      }
      res->noise_estimate = 0;
    }
  }

  // Estimate received pilot power
  if (isnormal(res->noise_estimate)) {
    res->snr = srslte_vec_avg_power_cf(q->pilot_recv_signal, nslots * nrefs_sym) / res->noise_estimate;
  } else {
    res->snr = NAN;
  }

  // Convert measurements in logarithm scale
  res->snr_db             = srslte_convert_power_to_dB(res->snr);
  res->noise_estimate_dbm = srslte_convert_power_to_dBm(res->noise_estimate);
}

int srslte_chest_ul_estimate_pusch(srslte_chest_ul_t*     q,
                                   srslte_ul_sf_cfg_t*    sf,
                                   srslte_pusch_cfg_t*    cfg,
                                   cf_t*                  input,
                                   srslte_chest_ul_res_t* res)
{
  if (!q->dmrs_signal_configured) {
    ERROR("Error must call srslte_chest_ul_set_cfg() before using the UL estimator\n");
    return SRSLTE_ERROR;
  }

  uint32_t nof_prb = cfg->grant.L_prb;

  if (!srslte_dft_precoding_valid_prb(nof_prb)) {
    ERROR("Error invalid nof_prb=%d\n", nof_prb);
    return SRSLTE_ERROR_INVALID_INPUTS;
  }

  int nrefs_sym = nof_prb * SRSLTE_NRE;
  int nrefs_sf  = nrefs_sym * SRSLTE_NOF_SLOTS_PER_SF;

  /* Get references from the input signal */
  srslte_refsignal_dmrs_pusch_get(&q->dmrs_signal, cfg, input, q->pilot_recv_signal); //copy input to pilot_recv_signal

  // Use the known DMRS signal to compute Least-squares estimates
  srslte_vec_prod_conj_ccc(q->pilot_recv_signal,
                           q->dmrs_pregen.r[cfg->grant.n_dmrs][sf->tti % SRSLTE_NOF_SF_X_FRAME][nof_prb],
                           q->pilot_estimates,
                           nrefs_sf);

  // Estimate
  chest_ul_estimate(q, SRSLTE_NOF_SLOTS_PER_SF, nrefs_sym, 1, cfg->meas_ta_en, true, cfg->grant.n_prb, res);

  return 0;
}

int srslte_chest_ul_estimate_pusch_with_tof(srslte_chest_ul_t*     q,
                                            srslte_ul_sf_cfg_t*    sf,
                                            srslte_pusch_cfg_t*    cfg,
                                            cf_t*                  input,
                                            srslte_chest_ul_res_t* res,
                                            float*                 tof,
                                            int                    save_to_file)
{
  if (!q->dmrs_signal_configured) {
    ERROR("Error must call srslte_chest_ul_set_cfg() before using the UL estimator\n");
    return SRSLTE_ERROR;
  }

  uint32_t nof_prb = cfg->grant.L_prb;

  if (!srslte_dft_precoding_valid_prb(nof_prb)) {
    ERROR("Error invalid nof_prb=%d\n", nof_prb);
    return SRSLTE_ERROR_INVALID_INPUTS;
  }

  int nrefs_sym = nof_prb * SRSLTE_NRE;
  int nrefs_sf  = nrefs_sym * SRSLTE_NOF_SLOTS_PER_SF;

  int num_extension = 3;

  /* Get references from the input signal */
  srslte_refsignal_dmrs_pusch_get(&q->dmrs_signal, cfg, input, q->pilot_recv_signal); //copy input to pilot_recv_signal


  { //do correlation and find the peak
    srslte_dft_plan_t ifft;
    srslte_dft_plan_c(&ifft, num_extension*nrefs_sf, SRSLTE_DFT_BACKWARD);

//    srslte_dft_plan_t fft;
//    srslte_dft_plan_c(&fft, nrefs_sf, SRSLTE_DFT_FORWARD);
//    *recv_buff = NULL, *dmrs_buff = NULL,
    cf_t *t_buff = NULL;
//    recv_buff = srslte_vec_cf_malloc(nrefs_sf);
//    dmrs_buff = srslte_vec_cf_malloc(nrefs_sf);
    t_buff    = srslte_vec_cf_malloc(num_extension*nrefs_sf);

//    srslte_dft_run_c(&fft, q->pilot_recv_signal, recv_buff);
//    srslte_dft_run_c(&fft, q->dmrs_pregen.r[cfg->grant.n_dmrs][sf->tti % SRSLTE_NOF_SF_X_FRAME][nof_prb], dmrs_buff);
    
    cf_t q_dmrs[num_extension*nrefs_sf];
    cf_t ref_dmrs[num_extension*nrefs_sf];
    
    cf_t *q_point = q->pilot_recv_signal;
    cf_t *ref_point = q->dmrs_pregen.r[cfg->grant.n_dmrs][sf->tti % SRSLTE_NOF_SF_X_FRAME][nof_prb];
    
    for(int i = 0; i < num_extension * nrefs_sf; i++) {
      if (i < nrefs_sf/2) {
        q_dmrs[i] = q_point[i];
        ref_dmrs[i] = ref_point[i];
      }
      else if (i >= (2*num_extension-1)*nrefs_sf/2) {
        q_dmrs[i] = q_point[i-(num_extension-1)*nrefs_sf];
        ref_dmrs[i] = ref_point[i-(num_extension-1)*nrefs_sf];
      }
      else {
        q_dmrs[i] = 0;
        ref_dmrs[i] = 0;
      }
    }
    // srslte_vec_prod_conj_ccc(q->pilot_recv_signal, q->dmrs_pregen.r[cfg->grant.n_dmrs][sf->tti % SRSLTE_NOF_SF_X_FRAME][nof_prb], t_buff, nrefs_sf);
    cf_t *q_pointer = q_dmrs;
    cf_t *ref_pointer = ref_dmrs;

    srslte_vec_prod_conj_ccc(q_pointer, ref_pointer, t_buff, num_extension*nrefs_sf);
    
    srslte_dft_run_c(&ifft, t_buff, t_buff);
    int max_ind = 0;
    float max_v = 0;
    for (int i = 0; i < num_extension*nrefs_sf; ++i) {
      if (cabsf(t_buff[i]) > max_v) {
        max_v = cabsf(t_buff[i]);
        max_ind = i;
      }
    }

    //printf("prn:%d, max ind: %d, nrefs_sf:%d\n",cfg->grant.L_prb, max_ind, nrefs_sf);
    //printf("ToF:%lf us\n", 2*66.67*(float)max_ind/(float)nrefs_sf);
    printf("???\n");
    printf("max_ind: %d\n", max_ind);
    printf("nrefs_sf: %d\n", nrefs_sf);
    printf("???\n");

    //write data to file
    if (save_to_file)
    {
      srslte_dft_plan_t ifft_0;
      srslte_dft_plan_c(&ifft_0, num_extension*nrefs_sf, SRSLTE_DFT_BACKWARD);
      cf_t q_dmrs_0[num_extension*nrefs_sf];
      cf_t ref_dmrs_0[num_extension*nrefs_sf];
      srslte_dft_run_c(&ifft_0, q_dmrs, q_dmrs_0);
      srslte_dft_run_c(&ifft_0, ref_dmrs, ref_dmrs_0);

      char file_name[80] = "\0";
      sprintf(file_name, "/home/dqs/data/data_%d.txt", max_ind);
      if (access(file_name, F_OK) != 0)
      {
        FILE *fptr = fopen(file_name, "w+");
        fprintf(fptr, "%d %d\n", num_extension*nrefs_sf, max_ind);
        for (int i = 0; i < num_extension*nrefs_sf; i++)
        {
          fprintf(fptr, "%f %f %f %f %f %f %f %f %f %f\n", crealf(q_pointer[i]), cimagf(q_pointer[i]),
                                                        crealf(q_dmrs_0[i]), cimagf(q_dmrs_0[i]),
                                                        crealf(ref_pointer[i]), cimagf(ref_pointer[i]),
                                                        crealf(ref_dmrs_0[i]), cimagf(ref_dmrs_0[i]),
                                                        crealf(t_buff[i]), cimagf(t_buff[i]));
        }
        fclose(fptr);
        printf("Data written to %s\n", file_name);
      }
    }

    if (tof != NULL) {
      if (max_ind < num_extension * nrefs_sf / 2) {
        *tof = 2*66.67*(float)max_ind/(float)(num_extension*nrefs_sf);
      }
      else {
        *tof = 2*66.67*(float)(max_ind-num_extension*nrefs_sf)/(float)(num_extension*nrefs_sf);
      }
    }
  }
  // Use the known DMRS signal to compute Least-squares estimates
  srslte_vec_prod_conj_ccc(q->pilot_recv_signal,
                           q->dmrs_pregen.r[cfg->grant.n_dmrs][sf->tti % SRSLTE_NOF_SF_X_FRAME][nof_prb],
                           q->pilot_estimates,
                           nrefs_sf);

  // Estimate
  chest_ul_estimate(q, SRSLTE_NOF_SLOTS_PER_SF, nrefs_sym, 1, cfg->meas_ta_en, true, cfg->grant.n_prb, res);

  return 0;
}


int srslte_chest_ul_estimate_pusch_debug(srslte_chest_ul_t*     q,
                                            srslte_ul_sf_cfg_t*    sf,
                                            srslte_pusch_cfg_t*    cfg,
                                            cf_t*                  input,
                                            srslte_chest_ul_res_t* res,
                                            float*                 tof,
                                            float*                 mean_val,
                                            float*                 mean_2_val,
                                            float*                 max_val,
                                            int                    save_to_file)
{
  if (!q->dmrs_signal_configured) {
    ERROR("Error must call srslte_chest_ul_set_cfg() before using the UL estimator\n");
    return SRSLTE_ERROR;
  }

  uint32_t nof_prb = cfg->grant.L_prb;

  if (!srslte_dft_precoding_valid_prb(nof_prb)) {
    ERROR("Error invalid nof_prb=%d\n", nof_prb);
    return SRSLTE_ERROR_INVALID_INPUTS;
  }

  int nrefs_sym = nof_prb * SRSLTE_NRE;
  int nrefs_sf  = nrefs_sym * SRSLTE_NOF_SLOTS_PER_SF;

  int num_extension = 3;

  /* Get references from the input signal */
  srslte_refsignal_dmrs_pusch_get(&q->dmrs_signal, cfg, input, q->pilot_recv_signal); //copy input to pilot_recv_signal


  { //do correlation and find the peak
    srslte_dft_plan_t ifft;
    srslte_dft_plan_c(&ifft, num_extension*nrefs_sf, SRSLTE_DFT_BACKWARD);

//    srslte_dft_plan_t fft;
//    srslte_dft_plan_c(&fft, nrefs_sf, SRSLTE_DFT_FORWARD);
//    *recv_buff = NULL, *dmrs_buff = NULL,
    cf_t *t_buff = NULL;
//    recv_buff = srslte_vec_cf_malloc(nrefs_sf);
//    dmrs_buff = srslte_vec_cf_malloc(nrefs_sf);
    t_buff    = srslte_vec_cf_malloc(num_extension*nrefs_sf);

//    srslte_dft_run_c(&fft, q->pilot_recv_signal, recv_buff);
//    srslte_dft_run_c(&fft, q->dmrs_pregen.r[cfg->grant.n_dmrs][sf->tti % SRSLTE_NOF_SF_X_FRAME][nof_prb], dmrs_buff);
    
    cf_t q_dmrs[num_extension*nrefs_sf];
    cf_t ref_dmrs[num_extension*nrefs_sf];
    
    cf_t *q_point = q->pilot_recv_signal;
    cf_t *ref_point = q->dmrs_pregen.r[cfg->grant.n_dmrs][sf->tti % SRSLTE_NOF_SF_X_FRAME][nof_prb];
    
    for(int i = 0; i < num_extension * nrefs_sf; i++) {
      if (i < nrefs_sf/2) {
        q_dmrs[i] = q_point[i];
        ref_dmrs[i] = ref_point[i];
      }
      else if (i >= (2*num_extension-1)*nrefs_sf/2) {
        q_dmrs[i] = q_point[i-(num_extension-1)*nrefs_sf];
        ref_dmrs[i] = ref_point[i-(num_extension-1)*nrefs_sf];
      }
      else {
        q_dmrs[i] = 0;
        ref_dmrs[i] = 0;
      }
    }
    // srslte_vec_prod_conj_ccc(q->pilot_recv_signal, q->dmrs_pregen.r[cfg->grant.n_dmrs][sf->tti % SRSLTE_NOF_SF_X_FRAME][nof_prb], t_buff, nrefs_sf);
    cf_t *q_pointer = q_dmrs;
    cf_t *ref_pointer = ref_dmrs;

    srslte_vec_prod_conj_ccc(q_pointer, ref_pointer, t_buff, num_extension*nrefs_sf);
    
    srslte_dft_run_c(&ifft, t_buff, t_buff);
    int max_ind = 0;
    float max_v = 0, mean_v = 0, mean_2_v = 0;
    for (int i = 0; i < num_extension*nrefs_sf; ++i) {
      if (cabsf(t_buff[i]) > max_v) {
        max_v = cabsf(t_buff[i]);
        max_ind = i;
      }
      mean_v += cabsf(t_buff[i]);
      mean_2_v += cabsf(t_buff[i])*cabsf(t_buff[i]);
    }
    mean_v /= (num_extension*nrefs_sf);
    mean_2_v /= (num_extension*nrefs_sf);

    //printf("prn:%d, max ind: %d, nrefs_sf:%d\n",cfg->grant.L_prb, max_ind, nrefs_sf);
    //printf("ToF:%lf us\n", 2*66.67*(float)max_ind/(float)nrefs_sf);
    printf("???\n");
    printf("max_ind: %d\n", max_ind);
    printf("nrefs_sf: %d\n", nrefs_sf);
    printf("???\n");

    //write data to file
    if (save_to_file)
    {
      srslte_dft_plan_t ifft_0;
      srslte_dft_plan_c(&ifft_0, num_extension*nrefs_sf, SRSLTE_DFT_BACKWARD);
      cf_t q_dmrs_0[num_extension*nrefs_sf];
      cf_t ref_dmrs_0[num_extension*nrefs_sf];
      srslte_dft_run_c(&ifft_0, q_dmrs, q_dmrs_0);
      srslte_dft_run_c(&ifft_0, ref_dmrs, ref_dmrs_0);

      char file_name[80] = "\0";
      sprintf(file_name, "/home/dqs/data/data_%d.txt", max_ind);
      if (access(file_name, F_OK) != 0)
      {
        FILE *fptr = fopen(file_name, "w+");
        fprintf(fptr, "%d %d\n", num_extension*nrefs_sf, max_ind);
        for (int i = 0; i < num_extension*nrefs_sf; i++)
        {
          fprintf(fptr, "%f %f %f %f %f %f %f %f %f %f\n", crealf(q_pointer[i]), cimagf(q_pointer[i]),
                                                        crealf(q_dmrs_0[i]), cimagf(q_dmrs_0[i]),
                                                        crealf(ref_pointer[i]), cimagf(ref_pointer[i]),
                                                        crealf(ref_dmrs_0[i]), cimagf(ref_dmrs_0[i]),
                                                        crealf(t_buff[i]), cimagf(t_buff[i]));
        }
        fclose(fptr);
        printf("Data written to %s\n", file_name);
      }
    }

    if (tof != NULL) {
      if (max_ind < num_extension * nrefs_sf / 2) {
        *tof = 2*66.67*(float)max_ind/(float)(num_extension*nrefs_sf);
      }
      else {
        *tof = 2*66.67*(float)(max_ind-num_extension*nrefs_sf)/(float)(num_extension*nrefs_sf);
      }
    }

    if (mean_val != NULL) {
      *mean_val = mean_v;
    }
    if (mean_2_val != NULL) {
      *mean_2_val = mean_2_v;
    }
    if (max_val != NULL) {
      *max_val = max_v;
    }

  }
  // Use the known DMRS signal to compute Least-squares estimates
  srslte_vec_prod_conj_ccc(q->pilot_recv_signal,
                           q->dmrs_pregen.r[cfg->grant.n_dmrs][sf->tti % SRSLTE_NOF_SF_X_FRAME][nof_prb],
                           q->pilot_estimates,
                           nrefs_sf);

  // Estimate
  chest_ul_estimate(q, SRSLTE_NOF_SLOTS_PER_SF, nrefs_sym, 1, cfg->meas_ta_en, true, cfg->grant.n_prb, res);

  return 0;
}


int srslte_chest_ul_estimate_pucch(srslte_chest_ul_t*     q,
                                   srslte_ul_sf_cfg_t*    sf,
                                   srslte_pucch_cfg_t*    cfg,
                                   cf_t*                  input,
                                   srslte_chest_ul_res_t* res)
{
  if (!q->dmrs_signal_configured) {
    ERROR("Error must call srslte_chest_ul_set_cfg() before using the UL estimator\n");
    return SRSLTE_ERROR;
  }

  int n_rs = srslte_refsignal_dmrs_N_rs(cfg->format, q->cell.cp);
  if (!n_rs) {
    ERROR("Error computing N_rs\n");
    return SRSLTE_ERROR;
  }
  int nrefs_sf = SRSLTE_NRE * n_rs * 2;

  /* Get references from the input signal */
  srslte_refsignal_dmrs_pucch_get(&q->dmrs_signal, cfg, input, q->pilot_recv_signal);

  /* Generate known pilots */
  if (cfg->format == SRSLTE_PUCCH_FORMAT_2A || cfg->format == SRSLTE_PUCCH_FORMAT_2B) {
    float max   = -1e9;
    int   i_max = 0;

    int m = 0;
    if (cfg->format == SRSLTE_PUCCH_FORMAT_2A) {
      m = 2;
    } else {
      m = 4;
    }

    for (int i = 0; i < m; i++) {
      cfg->pucch2_drs_bits[0] = i % 2;
      cfg->pucch2_drs_bits[1] = i / 2;
      srslte_refsignal_dmrs_pucch_gen(&q->dmrs_signal, sf, cfg, q->pilot_known_signal);
      srslte_vec_prod_conj_ccc(q->pilot_recv_signal, q->pilot_known_signal, q->pilot_estimates_tmp[i], nrefs_sf);
      float x = cabsf(srslte_vec_acc_cc(q->pilot_estimates_tmp[i], nrefs_sf));
      if (x >= max) {
        max   = x;
        i_max = i;
      }
    }
    memcpy(q->pilot_estimates, q->pilot_estimates_tmp[i_max], nrefs_sf * sizeof(cf_t));
    cfg->pucch2_drs_bits[0] = i_max % 2;
    cfg->pucch2_drs_bits[1] = i_max / 2;

  } else {
    srslte_refsignal_dmrs_pucch_gen(&q->dmrs_signal, sf, cfg, q->pilot_known_signal);
    /* Use the known DMRS signal to compute Least-squares estimates */
    srslte_vec_prod_conj_ccc(q->pilot_recv_signal, q->pilot_known_signal, q->pilot_estimates, nrefs_sf);
  }

  if (cfg->meas_ta_en && n_rs > 0) {
    float ta_err = 0.0;
    for (int ns = 0; ns < SRSLTE_NOF_SLOTS_PER_SF; ns++) {
      for (int i = 0; i < n_rs; i++) {
        ta_err += srslte_vec_estimate_frequency(&q->pilot_estimates[(i + ns * n_rs) * SRSLTE_NRE], SRSLTE_NRE) /
                  (float)(SRSLTE_NOF_SLOTS_PER_SF * n_rs);
      }
    }

    // Calculate actual time alignment error in micro-seconds
    if (isnormal(ta_err)) {
      ta_err /= 15e3f;                             // Convert from normalized frequency to seconds
      ta_err *= 1e6f;                              // Convert to micro-seconds
      ta_err     = roundf(ta_err * 10.0f) / 10.0f; // Round to one tenth of micro-second
      res->ta_us = ta_err;
    } else {
      res->ta_us = 0.0f;
    }
  }

  if (res->ce != NULL) {
    /* TODO: Currently averaging entire slot, performance good enough? */
    for (int ns = 0; ns < 2; ns++) {
      // Average all slot
      for (int i = 1; i < n_rs; i++) {
        srslte_vec_sum_ccc(&q->pilot_estimates[ns * n_rs * SRSLTE_NRE],
                           &q->pilot_estimates[(i + ns * n_rs) * SRSLTE_NRE],
                           &q->pilot_estimates[ns * n_rs * SRSLTE_NRE],
                           SRSLTE_NRE);
      }
      srslte_vec_sc_prod_ccc(&q->pilot_estimates[ns * n_rs * SRSLTE_NRE],
                             (float)1.0 / n_rs,
                             &q->pilot_estimates[ns * n_rs * SRSLTE_NRE],
                             SRSLTE_NRE);

      // Average in freq domain
      srslte_chest_average_pilots(&q->pilot_estimates[ns * n_rs * SRSLTE_NRE],
                                  &q->pilot_recv_signal[ns * n_rs * SRSLTE_NRE],
                                  q->smooth_filter,
                                  SRSLTE_NRE,
                                  1,
                                  q->smooth_filter_len);

      // Determine n_prb
      uint32_t n_prb = srslte_pucch_n_prb(&q->cell, cfg, ns);

      // copy estimates to slot
      for (int i = 0; i < SRSLTE_CP_NSYMB(q->cell.cp); i++) {
        memcpy(&res->ce[SRSLTE_RE_IDX(q->cell.nof_prb, i + ns * SRSLTE_CP_NSYMB(q->cell.cp), n_prb * SRSLTE_NRE)],
               &q->pilot_recv_signal[ns * n_rs * SRSLTE_NRE],
               sizeof(cf_t) * SRSLTE_NRE);
      }
    }
  }

  return 0;
}

int srslte_chest_ul_estimate_srs(srslte_chest_ul_t*                 q,
                                 srslte_ul_sf_cfg_t*                sf,
                                 srslte_refsignal_srs_cfg_t*        cfg,
                                 srslte_refsignal_dmrs_pusch_cfg_t* pusch_cfg,
                                 cf_t*                              input,
                                 srslte_chest_ul_res_t*             res)
{
  if (q == NULL || sf == NULL || cfg == NULL || pusch_cfg == NULL || input == NULL || res == NULL) {
    return SRSLTE_ERROR_INVALID_INPUTS;
  }

  // Extract parameters
  uint32_t n_srs_re = srslte_refsignal_srs_M_sc(&q->dmrs_signal, cfg);

  // Extract Sounding Reference Signal
  if (srslte_refsignal_srs_get(&q->dmrs_signal, cfg, sf->tti, q->pilot_recv_signal, input) != SRSLTE_SUCCESS) {
    return SRSLTE_ERROR;
  }

  // Get Known pilots
  cf_t* known_pilots = q->pilot_known_signal;
  if (q->srs_signal_configured) {
    known_pilots = q->srs_pregen.r[sf->tti % SRSLTE_NOF_SF_X_FRAME];
  } else {
    srslte_refsignal_srs_gen(&q->dmrs_signal, cfg, pusch_cfg, sf->tti % SRSLTE_NOF_SF_X_FRAME, known_pilots);
  }

  // Compute least squares estimates
  srslte_vec_prod_conj_ccc(q->pilot_recv_signal, known_pilots, q->pilot_estimates, n_srs_re);

  // Estimate
  uint32_t n_prb[2] = {};
  chest_ul_estimate(q, 1, n_srs_re, 1, true, false, n_prb, res);

  return SRSLTE_SUCCESS;
}