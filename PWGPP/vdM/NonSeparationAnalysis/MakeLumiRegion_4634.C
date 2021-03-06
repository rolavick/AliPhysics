// -*- C++ -*-
// $Id$

void MakeLumiRegion_4634() {
  gROOT->LoadMacro("AliLuminousRegionFit.cxx+");
  gROOT->LoadMacro("Util.C");

  const ScanData sd[] = {
    { 4634,
      11,
      "root/4634/AnalysisResults_VdM_244369.root",
      "txt/4634/SepVStime.out",
      "Scan1",
      "Scan1X", 0, 1447981655, 1447982680, 0,
      "Scan1Y", 1, 1447982863, 1447983886, 0                             
    },
    
    { 4634,
      11,
      "root/4634/AnalysisResults_VdM_244369.root",
      "txt/4634/SepVStime.out",
      "Scan2",
      "Scan2X", 0, 1447984119, 1447985146, 0,
      "Scan2Y", 1, 1447985327, 1447986351, 0
    },
    
    { 4634,
      11,
      "root/4634/AnalysisResults_VdM_244375.root",
      "txt/4634/SepVStime.out",
      "ScanOffset",
      "ScanOffsetX", 0, 1447988540, 1447989542, -457.5,
      "ScanOffsetY", 1, 1447989790, 1447990716, -457.5
    }
  };

  const Int_t n = sizeof(sd)/sizeof(ScanData);

  const Int_t bcs[] = {
      -1, // -B
     344, // -I1 BCM5   344H1L3219H
     464, // -I2 BCM6   464H1L3099H
     827, // -I3 BCM7   827H1L2736H
    1187, // -I4 BCM8  1187H1L2376H
    1558, // -I5 BCM9  1558H1L2005H
    1678, // -I6 BCM10 1678H1L1885H
    3177, // -I7 BCM11 3177H1L386H
    3297  // -I8 BCM12 3297H1L266H
  };

  for (Int_t i=0; i<n; ++i) {
    CheckCopyFile(sd[i].vtxFileName);
    AliLuminousRegionFit f(sd[i].fillNumber,
			   sd[i].minNumberOfTracks,
			   sd[i].vtxFileName,
			   sd[i].sepFileName);
    for (Int_t j=0, m=sizeof(bcs)/sizeof(Int_t); j<m; ++j) {
      f.DoFit(sd[i].scanName1, sd[i].t1, sd[i].t2, sd[i].scanType1, sd[i].offset1, bcs[j]);
      f.DoFit(sd[i].scanName2, sd[i].t3, sd[i].t4, sd[i].scanType2, sd[i].offset2, bcs[j]);
    }
  }
}
