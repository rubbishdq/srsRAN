/**
 *
 * \section COPYRIGHT
 *
 * Copyright 2013-2020 Software Radio Systems Limited
 *
 * By using this file, you agree to the terms and conditions set
 * forth in the LICENSE file which can be found at the top level of
 * the distribution.
 *
 */

#ifndef SRSLTE_MAC_CONTROLLER_H
#define SRSLTE_MAC_CONTROLLER_H

#include "rrc_bearer_cfg.h"
#include "rrc_cell_cfg.h"
#include "srslte/interfaces/rrc_interface_types.h"
#include "srslte/interfaces/sched_interface.h"
#include <bitset>

namespace srsenb {

class mac_interface_rrc;

class mac_controller
{
  using ue_cfg_t = sched_interface::ue_cfg_t;

public:
  mac_controller(uint16_t                    rnti_,
                 const ue_cell_ded_list&     ue_cell_list_,
                 const bearer_cfg_handler&   bearer_list_,
                 const rrc_cfg_t&            rrc_cfg_,
                 mac_interface_rrc*          mac_,
                 const enb_cell_common_list& cell_common_list,
                 const ue_cfg_t&             sched_ue_cfg);

  // Handling of Msg4
  int  handle_con_setup(const asn1::rrc::rrc_conn_setup_r8_ies_s& conn_setup);
  int  handle_con_reest(const asn1::rrc::rrc_conn_reest_r8_ies_s& conn_reest);
  void handle_con_reject();
  int  handle_crnti_ce(uint32_t temp_crnti);

  void handle_con_setup_complete();
  void handle_con_reest_complete();

  void handle_con_reconf(const asn1::rrc::rrc_conn_recfg_r8_ies_s& conn_recfg,
                         const srslte::rrc_ue_capabilities_t&      uecaps);
  void handle_con_reconf_complete();

  void handle_target_enb_ho_cmd(const asn1::rrc::rrc_conn_recfg_r8_ies_s& conn_recfg,
                                const srslte::rrc_ue_capabilities_t&      uecaps);
  void handle_intraenb_ho_cmd(const asn1::rrc::rrc_conn_recfg_r8_ies_s& conn_recfg,
                              const srslte::rrc_ue_capabilities_t&      uecaps);
  void handle_ho_prep(const asn1::rrc::ho_prep_info_r8_ies_s& ho_prep);

  void handle_max_retx();

  const ue_cfg_t& get_ue_sched_cfg() const { return current_sched_ue_cfg; }
  bool            is_crnti_set() const { return crnti_set; }

  void set_scell_activation(const std::bitset<SRSLTE_MAX_CARRIERS>& scell_mask);
  void set_drb_activation(bool active);

  enum proc_stage_t : int8_t { config_tx, config_complete, other };
  void update_mac(proc_stage_t stage);

private:
  int  apply_basic_conn_cfg(const asn1::rrc::rr_cfg_ded_s& rr_cfg);
  void apply_current_bearers_cfg();

  srslog::basic_logger&       logger;
  uint16_t                    rnti;
  const ue_cell_ded_list&     ue_cell_list;
  const bearer_cfg_handler&   bearer_list;
  const rrc_cfg_t*            rrc_cfg = nullptr;
  mac_interface_rrc*          mac     = nullptr;
  const enb_cell_common_list& cell_common_list;

  /// UE configuration currently present at the MAC, including any transient disabling of bearers/scells
  ue_cfg_t current_sched_ue_cfg = {};
  /// UE configuration once the RRC config procedure (e.g. Reconfiguration) is complete
  ue_cfg_t next_sched_ue_cfg = {};
  bool     crnti_set         = false;
};

} // namespace srsenb

#endif // SRSLTE_MAC_CONTROLLER_H
