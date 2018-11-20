#!/bin/bash
mkdir etc
mkdir plots
wget -O etc/rtcoef_metop_2_iasi_so2.H5 https://www.nwpsaf.eu/downloads/rtcoef_rttov12/rttov9pred101L/rtcoef_metop_2_iasi_so2.H5
wget -O etc/rtcoef_eos_2_airs_so2.H5 https://www.nwpsaf.eu/downloads/rtcoef_rttov12/rttov9pred101L/rtcoef_eos_2_airs_so2.H5
wget -O etc/rtcoef_jpss_0_cris_so2.H5 https://www.nwpsaf.eu/downloads/rtcoef_rttov12/rttov9pred101L/rtcoef_jpss_0_cris_so2.H5
wget -O etc/rtcoef_jpss_0_cris-fsr_so2.H5 https://www.nwpsaf.eu/downloads/rtcoef_rttov12/rttov9pred101L/rtcoef_jpss_0_cris-fsr_so2.H5

dropbox_uploader.sh download d07bin_oc_day.f89 etc/d07bin_oc_day.f89
dropbox_uploader.sh download airs_v5p9.f89 etc/airs_v5p9.f89
dropbox_uploader.sh download cris_npp_fsr.f89  etc/cris_npp_fsr.f89
dropbox_uploader.sh download cris_npp_nsr.f89 etc/cris_npp_nsr.f89
dropbox_uploader.sh download active_channels.tbl etc/active_channels.tbl
