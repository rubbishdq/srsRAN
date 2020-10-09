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

#ifndef SRSLTE_MAC_CONTROLLER_H
#define SRSLTE_MAC_CONTROLLER_H

#include "rrc_ue.h"
#include <bitset>

namespace srsenb {

class rrc::ue::mac_controller
{
public:
  mac_controller(rrc::ue* rrc_ue, const sched_interface::ue_cfg_t& sched_ue_cfg);

  // Handling of Msg4
  int  handle_con_setup(const asn1::rrc::rrc_conn_setup_r8_ies_s& conn_setup);
  int  handle_con_reest(const asn1::rrc::rrc_conn_reest_r8_ies_s& conn_reest);
  void handle_con_reject();
  int  handle_crnti_ce(uint32_t temp_crnti);

  void handle_con_setup_complete();
  void handle_con_reest_complete();

  void handle_con_reconf(const asn1::rrc::rrc_conn_recfg_r8_ies_s& conn_recfg);
  void handle_con_reconf_complete();

  void handle_ho_prep(const asn1::rrc::ho_prep_info_r8_ies_s&   ho_prep,
                      const asn1::rrc::rrc_conn_recfg_r8_ies_s& conn_recfg);
  void handle_ho_prep_complete();

  const sched_interface::ue_cfg_t& get_ue_sched_cfg() const { return current_sched_ue_cfg; }

  void set_scell_activation(const std::bitset<SRSLTE_MAX_CARRIERS>& scell_mask);

  enum proc_stage_t : int8_t { config_tx, config_complete, other };
  void update_mac(proc_stage_t stage);

private:
  void handle_con_reconf_with_mobility();
  int  apply_basic_conn_cfg(const asn1::rrc::rr_cfg_ded_s& rr_cfg);
  void apply_current_bearers_cfg();
  void apply_phy_cfg_updates_common(const asn1::rrc::phys_cfg_ded_s& phy_cfg);
  void apply_scell_cfg_updates(uint32_t ue_cc_idx);

  srslte::log_ref           log_h;
  rrc_cfg_t*                rrc_cfg              = nullptr;
  rrc::ue*                  rrc_ue               = nullptr;
  mac_interface_rrc*        mac                  = nullptr;
  sched_interface::ue_cfg_t current_sched_ue_cfg = {};
  bool                      crnti_set            = false;

  // pending changes
  std::unique_ptr<asn1::rrc::scell_to_add_mod_list_r10_l> pending_scells_cfg;
};

} // namespace srsenb

#endif // SRSLTE_MAC_CONTROLLER_H
