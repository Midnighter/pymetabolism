Maximize
  + Biomass_Ecoli_core_N__w_GAM_
Subject To
 glc_D_e: EX_glc_e__Rev - EX_glc_e_ - GLCpts = 0
 glc_D_b: - EX_glc_e__Rev + EX_glc_e_ + glc_D_b_Source - glc_D_b_Drain = 0
 pep_c: - GLCpts - PYK - PPC - ENO_Rev + ENO + PPS
   - 0.5191 Biomass_Ecoli_core_N__w_GAM_ + PPCK = 0
 pyr_c: GLCpts - LDH_D_Rev + LDH_D + ME1 + ME2 - PFL + PYK - PPS
   - 2.8328 Biomass_Ecoli_core_N__w_GAM_ - PDH - PYRt2r_Rev + PYRt2r = 0
 g6p_c: GLCpts + PGI_Rev - PGI + G6PDH2r_Rev - G6PDH2r
   - 0.205 Biomass_Ecoli_core_N__w_GAM_ = 0
 co2_e: EX_co2_e__Rev - EX_co2_e_ + CO2t_Rev - CO2t = 0
 co2_b: - EX_co2_e__Rev + EX_co2_e_ + co2_b_Source - co2_b_Drain = 0
 for_e: - EX_for_e_ + FORt_Rev - FORt = 0
 for_b: EX_for_e_ + for_b_Source - for_b_Drain = 0
 r5p_c: RPI_Rev - RPI + TKT1_Rev - TKT1
   - 0.8977 Biomass_Ecoli_core_N__w_GAM_ = 0
 ru5p_D_c: - RPI_Rev + RPI + GND + RPE_Rev - RPE = 0
 xu5p_D_c: TKT1_Rev - TKT1 + TKT2_Rev - TKT2 - RPE_Rev + RPE = 0
 g3p_c: - TKT1_Rev + TKT1 - TPI_Rev + TPI + GAPD_Rev - GAPD - FBA_Rev + FBA
   - 0.129 Biomass_Ecoli_core_N__w_GAM_ - TKT2_Rev + TKT2 + TALA_Rev - TALA
   = 0
 s7p_c: - TKT1_Rev + TKT1 + TALA_Rev - TALA = 0
 3pg_c: PGK_Rev - PGK - PGM_Rev + PGM - 1.496 Biomass_Ecoli_core_N__w_GAM_
   = 0
 atp_c: PGK_Rev - PGK + SUCOAS_Rev - SUCOAS + ADK1_Rev - ADK1 + PYK
   + ACKr_Rev - ACKr - PFK - PPS - ATPS4r_Rev + ATPS4r - ATPM
   - 55.703 Biomass_Ecoli_core_N__w_GAM_ - PPCK = 0
 adp_c: - PGK_Rev + PGK - SUCOAS_Rev + SUCOAS - 2 ADK1_Rev + 2 ADK1 - PYK
   - ACKr_Rev + ACKr + PFK + ATPS4r_Rev - ATPS4r + ATPM
   + 55.703 Biomass_Ecoli_core_N__w_GAM_ + PPCK = 0
 13dpg_c: - PGK_Rev + PGK - GAPD_Rev + GAPD = 0
 h2o_e: EX_h2o_e__Rev - EX_h2o_e_ + H2Ot_Rev - H2Ot = 0
 h2o_b: - EX_h2o_e__Rev + EX_h2o_e_ + h2o_b_Source - h2o_b_Drain = 0
 nadp_c: - GND + G6PDH2r_Rev - G6PDH2r - ME2 - THD2 + NADTRHD + ICDHyr_Rev
   - ICDHyr + 18.225 Biomass_Ecoli_core_N__w_GAM_ = 0
 6pgc_c: - GND + PGL = 0
 nadph_c: GND - G6PDH2r_Rev + G6PDH2r + ME2 + THD2 - NADTRHD - ICDHyr_Rev
   + ICDHyr - 18.225 Biomass_Ecoli_core_N__w_GAM_ = 0
 co2_c: GND - CO2t_Rev + CO2t + ME1 + ME2 + AKGDH - PPC - ICDHyr_Rev
   + ICDHyr + PDH + PPCK = 0
 fum_e: - EX_fum_e_ - FUMt2_2 = 0
 fum_b: EX_fum_e_ + fum_b_Source - fum_b_Drain = 0
 succ_c: SUCOAS_Rev - SUCOAS + SUCCt2_2 - SUCCt2b - SUCD1i + FRD + ICL = 0
 coa_c: SUCOAS_Rev - SUCOAS + MALS - PTAr_Rev + PTAr - ADHEr_Rev + ADHEr
   - PFL - AKGDH + 3.7478 Biomass_Ecoli_core_N__w_GAM_ + CS - PDH = 0
 pi_c: - SUCOAS_Rev + SUCOAS + PTAr_Rev - PTAr + PIt_Rev - PIt + PPC + FBP
   + PPS + GAPD_Rev - GAPD + ATPS4r_Rev - ATPS4r + ATPM
   + 55.703 Biomass_Ecoli_core_N__w_GAM_ = 0
 succoa_c: - SUCOAS_Rev + SUCOAS + AKGDH = 0
 lac_D_c: LDH_D_Rev - LDH_D - D_LACt2_Rev + D_LACt2 = 0
 nad_c: LDH_D_Rev - LDH_D - 2 ADHEr_Rev + 2 ADHEr + MDH_Rev - MDH - ME1
   + THD2 - AKGDH + NADH11 - NADTRHD + GAPD_Rev - GAPD
   - 3.547 Biomass_Ecoli_core_N__w_GAM_ - PDH = 0
 h_c: - LDH_D_Rev + LDH_D + 2 SUCCt2_2 + MALS + PGL - ACt2r_Rev + ACt2r
   - G6PDH2r_Rev + G6PDH2r - D_LACt2_Rev + D_LACt2 + 2 ADHEr_Rev - 2 ADHEr
   - MDH_Rev + MDH - SUCCt2b + 2 THD2 - AKGt2r_Rev + AKGt2r - 3 NADH11
   - PYK + PPC + PFK + 2 PPS - GAPD_Rev + GAPD - 3 ATPS4r_Rev + 3 ATPS4r
   + 2 FUMt2_2 + ATPM + 41.025 Biomass_Ecoli_core_N__w_GAM_ + CS - 2 CYTBD
   - ETOHt2r_Rev + ETOHt2r - PYRt2r_Rev + PYRt2r = 0
 nadh_c: - LDH_D_Rev + LDH_D + 2 ADHEr_Rev - 2 ADHEr - MDH_Rev + MDH + ME1
   - THD2 + AKGDH - NADH11 + NADTRHD - GAPD_Rev + GAPD
   + 3.547 Biomass_Ecoli_core_N__w_GAM_ + PDH = 0
 h_e: EX_h_e__Rev - EX_h_e_ - 2 SUCCt2_2 + ACt2r_Rev - ACt2r + D_LACt2_Rev
   - D_LACt2 + SUCCt2b - 2 THD2 + AKGt2r_Rev - AKGt2r + 2 NADH11
   + 4 ATPS4r_Rev - 4 ATPS4r - 2 FUMt2_2 + 2 CYTBD + ETOHt2r_Rev - ETOHt2r
   + PYRt2r_Rev - PYRt2r = 0
 h_b: - EX_h_e__Rev + EX_h_e_ + h_b_Source - h_b_Drain = 0
 etoh_e: - EX_etoh_e_ + ETOHt2r_Rev - ETOHt2r = 0
 etoh_b: EX_etoh_e_ + etoh_b_Source - etoh_b_Drain = 0
 h2o_c: - H2Ot_Rev + H2Ot - MALS - PGL - PPC - FBP - ENO_Rev + ENO - PPS
   - ATPS4r_Rev + ATPS4r - ATPM - 55.703 Biomass_Ecoli_core_N__w_GAM_
   + FUM_Rev - FUM - CS + CYTBD = 0
 o2_e: O2t_Rev - O2t + EX_o2_e__Rev - EX_o2_e_ = 0
 o2_c: - O2t_Rev + O2t - 0.5 CYTBD = 0
 succ_e: - SUCCt2_2 + SUCCt2b - EX_succ_e_ = 0
 glx_c: - MALS + ICL = 0
 accoa_c: - MALS + PTAr_Rev - PTAr + ADHEr_Rev - ADHEr + PFL
   - 3.7478 Biomass_Ecoli_core_N__w_GAM_ - CS + PDH = 0
 mal_L_c: MALS + MDH_Rev - MDH - ME1 - ME2 - FUM_Rev + FUM = 0
 actp_c: - PTAr_Rev + PTAr - ACKr_Rev + ACKr = 0
 q8_c: SUCD4_Rev - SUCD4 - NADH11 + CYTBD = 0
 fadh2_c: SUCD4_Rev - SUCD4 + SUCD1i - FRD = 0
 fad_c: - SUCD4_Rev + SUCD4 - SUCD1i + FRD = 0
 q8h2_c: - SUCD4_Rev + SUCD4 + NADH11 - CYTBD = 0
 pi_e: EX_pi_e__Rev - EX_pi_e_ - PIt_Rev + PIt = 0
 pi_b: - EX_pi_e__Rev + EX_pi_e_ + pi_b_Source - pi_b_Drain = 0
 6pgl_c: - PGL - G6PDH2r_Rev + G6PDH2r = 0
 dhap_c: TPI_Rev - TPI - FBA_Rev + FBA = 0
 o2_b: - EX_o2_e__Rev + EX_o2_e_ + o2_b_Source - o2_b_Drain = 0
 ac_e: ACt2r_Rev - ACt2r - EX_ac_e_ = 0
 ac_c: - ACt2r_Rev + ACt2r + ACKr_Rev - ACKr = 0
 2pg_c: PGM_Rev - PGM + ENO_Rev - ENO = 0
 pyr_e: - EX_pyr_e_ + PYRt2r_Rev - PYRt2r = 0
 pyr_b: EX_pyr_e_ + pyr_b_Source - pyr_b_Drain = 0
 f6p_c: - PGI_Rev + PGI + FBP - PFK - 0.0709 Biomass_Ecoli_core_N__w_GAM_
   - TKT2_Rev + TKT2 - TALA_Rev + TALA = 0
 lac_D_e: D_LACt2_Rev - D_LACt2 - EX_lac_D_e_ = 0
 etoh_c: - ADHEr_Rev + ADHEr - ETOHt2r_Rev + ETOHt2r = 0
 oaa_c: - MDH_Rev + MDH + PPC - 1.7867 Biomass_Ecoli_core_N__w_GAM_ - CS
   - PPCK = 0
 lac_D_b: EX_lac_D_e_ + lac_D_b_Source - lac_D_b_Drain = 0
 fum_c: SUCD1i - FRD + FUMt2_2 + FUM_Rev - FUM = 0
 amp_c: ADK1_Rev - ADK1 + PPS = 0
 akg_e: AKGt2r_Rev - AKGt2r - EX_akg_e_ = 0
 akg_c: - AKGt2r_Rev + AKGt2r - AKGDH - ICDHyr_Rev + ICDHyr
   - 1.0789 Biomass_Ecoli_core_N__w_GAM_ = 0
 for_c: PFL - FORt_Rev + FORt = 0
 fdp_c: - FBP + PFK + FBA_Rev - FBA = 0
 icit_c: ICDHyr_Rev - ICDHyr - ACONT_Rev + ACONT - ICL = 0
 akg_b: EX_akg_e_ + akg_b_Source - akg_b_Drain = 0
 cit_c: ACONT_Rev - ACONT + CS = 0
 e4p_c: - 0.361 Biomass_Ecoli_core_N__w_GAM_ + TKT2_Rev - TKT2 - TALA_Rev
   + TALA = 0
 ac_b: EX_ac_e_ + ac_b_Source - ac_b_Drain = 0
 succ_b: EX_succ_e_ + succ_b_Source - succ_b_Drain = 0
Bounds
 EX_glc_e__Rev <= 10
 EX_glc_e_ <= 999999
 GLCpts <= 999999
 EX_co2_e__Rev <= 999999
 EX_co2_e_ <= 999999
 EX_for_e_ <= 999999
 RPI_Rev <= 999999
 RPI <= 999999
 TKT1_Rev <= 999999
 TKT1 <= 999999
 PGK_Rev <= 999999
 PGK <= 999999
 EX_h2o_e__Rev <= 999999
 EX_h2o_e_ <= 999999
 GND <= 999999
 EX_fum_e_ <= 999999
 SUCOAS_Rev <= 999999
 SUCOAS <= 999999
 LDH_D_Rev <= 999999
 LDH_D <= 999999
 EX_h_e__Rev <= 999999
 EX_h_e_ <= 999999
 EX_etoh_e_ <= 999999
 H2Ot_Rev <= 999999
 H2Ot <= 999999
 O2t_Rev <= 999999
 O2t <= 999999
 SUCCt2_2 <= 999999
 MALS <= 999999
 PTAr_Rev <= 999999
 PTAr <= 999999
 SUCD4_Rev <= 999999
 SUCD4 <= 999999
 EX_pi_e__Rev <= 999999
 EX_pi_e_ <= 999999
 PGL <= 999999
 TPI_Rev <= 999999
 TPI <= 999999
 EX_o2_e__Rev <= 999999
 EX_o2_e_ <= 999999
 ACt2r_Rev <= 999999
 ACt2r <= 999999
 CO2t_Rev <= 999999
 CO2t <= 999999
 PGM_Rev <= 999999
 PGM <= 999999
 EX_pyr_e_ <= 999999
 PGI_Rev <= 999999
 PGI <= 999999
 G6PDH2r_Rev <= 999999
 G6PDH2r <= 999999
 D_LACt2_Rev <= 999999
 D_LACt2 <= 999999
 ADHEr_Rev <= 999999
 ADHEr <= 999999
 MDH_Rev <= 999999
 MDH <= 999999
 EX_lac_D_e_ <= 999999
 SUCCt2b <= 999999
 ME1 <= 999999
 PIt_Rev <= 999999
 PIt <= 999999
 SUCD1i <= 999999
 ME2 <= 999999
 FRD <= 999999
 THD2 <= 999999
 ADK1_Rev <= 999999
 ADK1 <= 999999
 AKGt2r_Rev <= 999999
 AKGt2r <= 999999
 PFL <= 999999
 AKGDH <= 999999
 NADH11 <= 999999
 NADTRHD <= 999999
 PYK <= 999999
 PPC <= 999999
 FBP <= 999999
 ACKr_Rev <= 999999
 ACKr <= 999999
 ENO_Rev <= 999999
 ENO <= 999999
 PFK <= 999999
 FORt_Rev <= 999999
 FORt <= 999999
 ICDHyr_Rev <= 999999
 ICDHyr <= 999999
 EX_akg_e_ <= 999999
 PPS <= 999999
 ACONT_Rev <= 999999
 ACONT <= 999999
 GAPD_Rev <= 999999
 GAPD <= 999999
 ATPS4r_Rev <= 999999
 ATPS4r <= 999999
 FUMt2_2 <= 999999
 FBA_Rev <= 999999
 FBA <= 999999
 ATPM = 7.6
 Biomass_Ecoli_core_N__w_GAM_ <= 999999
 FUM_Rev <= 999999
 FUM <= 999999
 CS <= 999999
 PDH <= 999999
 ICL <= 999999
 CYTBD <= 999999
 PPCK <= 999999
 EX_ac_e_ <= 999999
 TKT2_Rev <= 999999
 TKT2 <= 999999
 EX_succ_e_ <= 999999
 RPE_Rev <= 999999
 RPE <= 999999
 ETOHt2r_Rev <= 999999
 ETOHt2r <= 999999
 TALA_Rev <= 999999
 TALA <= 999999
 PYRt2r_Rev <= 999999
 PYRt2r <= 999999
 ac_b_Source <= 1000
 ac_b_Drain <= 1000
 glc_D_b_Source <= 1000
 glc_D_b_Drain <= 1000
 pi_b_Source <= 1000
 pi_b_Drain <= 1000
 o2_b_Source <= 1000
 o2_b_Drain <= 1000
 etoh_b_Source <= 1000
 etoh_b_Drain <= 1000
 h2o_b_Source <= 1000
 h2o_b_Drain <= 1000
 pyr_b_Source <= 1000
 pyr_b_Drain <= 1000
 h_b_Source <= 1000
 h_b_Drain <= 1000
 fum_b_Source <= 1000
 fum_b_Drain <= 1000
 akg_b_Source <= 1000
 akg_b_Drain <= 1000
 lac_D_b_Source <= 1000
 lac_D_b_Drain <= 1000
 for_b_Source <= 1000
 for_b_Drain <= 1000
 co2_b_Source <= 1000
 co2_b_Drain <= 1000
 succ_b_Source <= 1000
 succ_b_Drain <= 1000
End
