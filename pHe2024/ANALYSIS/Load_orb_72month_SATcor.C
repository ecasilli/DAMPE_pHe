{
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
//
//
TChain *skim= new TChain("newtree");

// =============================
//---------Orbital data------
//---------Orbital data------
//Jan 2016
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201601/skim_data_201601_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201601/skim_data_201601_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201601/skim_data_201601_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Jan_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Jan_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201601/skim_data_201601_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201601/skim_data_201601_part2_002_010.root");
//Feb 2016
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201602/skim_data_201602_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201602/skim_data_201602_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201602/skim_data_201602_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Feb_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Feb_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201602/skim_data_201602_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201602/skim_data_201602_part2_002_010.root");
// Mar 2016
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201603/skim_data_201603_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201603/skim_data_201603_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201603/skim_data_201603_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Mar_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Mar_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201603/skim_data_201603_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201603/skim_data_201603_part2_002_010.root");
//

// Apr 2016
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201604/skim_data_201604_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201604/skim_data_201604_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201604/skim_data_201604_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Apr_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Apr_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201604/skim_data_201604_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201604/skim_data_201604_part2_002_010.root");
//
// May 2016
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201605/skim_data_201605_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201605/skim_data_201605_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201605/skim_data_201605_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016May_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016May_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201605/skim_data_201605_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201605/skim_data_201605_part2_002_010.root");
//
// Jun 2016
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201606/skim_data_201606_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201606/skim_data_201606_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201606/skim_data_201606_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Jun_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Jun_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201606/skim_data_201606_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201606/skim_data_201606_part2_002_010.root");
//
// Jul 2016
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201607/skim_data_201607_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201607/skim_data_201607_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201607/skim_data_201607_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Jul_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Jul_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201607/skim_data_201607_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201607/skim_data_201607_part2_002_010.root");
//
// Aug 2016
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201608/skim_data_201608_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201608/skim_data_201608_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201608/skim_data_201608_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Aug_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Aug_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201608/skim_data_201608_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201608/skim_data_201608_part2_002_010.root");
//
// Sep 2016
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201609/skim_data_201609_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201609/skim_data_201609_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201609/skim_data_201609_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Sep_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Sep_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201609/skim_data_201609_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201609/skim_data_201609_part2_002_010.root");
//
// Oct 2016
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201610/skim_data_201610_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201610/skim_data_201610_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201610/skim_data_201610_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Oct_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Oct_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201610/skim_data_201610_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201610/skim_data_201610_part2_002_010.root");
//
// Nov 2016
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201611/skim_data_201611_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201611/skim_data_201611_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201611/skim_data_201611_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Nov_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Nov_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201611/skim_data_201611_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201611/skim_data_201611_part2_002_010.root");
//
// Dec 2016
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201612/skim_data_201612_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201612/skim_data_201612_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201612/skim_data_201612_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Dec_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2016Dec_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201612/skim_data_201612_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201612/skim_data_201612_part2_002_010.root");
//
// --------------------
//
// Jan 2017
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201701/skim_data_201701_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201701/skim_data_201701_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201701/skim_data_201701_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Jan_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Jan_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201701/skim_data_201701_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201701/skim_data_201701_part2_002_010.root");
//
// Feb 2017
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201702/skim_data_201702_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201702/skim_data_201702_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201702/skim_data_201702_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Feb_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Feb_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201702/skim_data_201702_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201702/skim_data_201702_part2_002_010.root");
//
// Mar 2017
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201703/skim_data_201703_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201703/skim_data_201703_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201703/skim_data_201703_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Mar_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Mar_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201703/skim_data_201703_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201703/skim_data_201703_part2_002_010.root");
//
// Apr 2017
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201704/skim_data_201704_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201704/skim_data_201704_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201704/skim_data_201704_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Apr_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Apr_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201704/skim_data_201704_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201704/skim_data_201704_part2_002_010.root");
//
// May 2017
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201705/skim_data_201705_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201705/skim_data_201705_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201705/skim_data_201705_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017May_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017May_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201705/skim_data_201705_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201705/skim_data_201705_part2_002_010.root");
//
// Jun 2017
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201706/skim_data_201706_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201706/skim_data_201706_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201706/skim_data_201706_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Jun_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Jun_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201706/skim_data_201706_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201706/skim_data_201706_part2_002_010.root");
//
// Jul 2017
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201707/skim_data_201707_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201707/skim_data_201707_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201707/skim_data_201707_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Jul_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Jul_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201707/skim_data_201707_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201707/skim_data_201707_part2_002_010.root");
//
// Aug 2017
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201708/skim_data_201708_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201708/skim_data_201708_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201708/skim_data_201708_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Aug_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Aug_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201708/skim_data_201708_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201708/skim_data_201708_part2_002_010.root");
//
// Sep 2017
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201709/skim_data_201709_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201709/skim_data_201709_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201709/skim_data_201709_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Sep_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Sep_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201709/skim_data_201709_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201709/skim_data_201709_part2_002_010.root");
//
// Oct 2017
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201710/skim_data_201710_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201710/skim_data_201710_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201710/skim_data_201710_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Oct_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Oct_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201710/skim_data_201710_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201710/skim_data_201710_part2_002_010.root");

// Nov 2017
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201711/skim_data_201711_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201711/skim_data_201711_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201711/skim_data_201711_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Nov_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Nov_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201711/skim_data_201711_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201711/skim_data_201711_part2_002_010.root");


// Dec 2017
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201712/skim_data_201712_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201712/skim_data_201712_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201712/skim_data_201712_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Dec_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2017Dec_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201712/skim_data_201712_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201712/skim_data_201712_part2_002_010.root");

//Jan 2018
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201801/skim_data_201801_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201801/skim_data_201801_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201801/skim_data_201801_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Jan_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Jan_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201801/skim_data_201801_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201801/skim_data_201801_part2_002_010.root");

//Feb 2018
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201802/skim_data_201802_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201802/skim_data_201802_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201802/skim_data_201802_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Feb_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Feb_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201802/skim_data_201802_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201802/skim_data_201802_part2_002_010.root");

//Mar 2018
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201803/skim_data_201803_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201803/skim_data_201803_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201803/skim_data_201803_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Mar_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Mar_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201803/skim_data_201803_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201803/skim_data_201803_part2_002_010.root");


//Apr 2018
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201804/skim_data_201804_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201804/skim_data_201804_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201804/skim_data_201804_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Apr_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Apr_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201804/skim_data_201804_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201804/skim_data_201804_part2_002_010.root");

//May 2018
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201805/skim_data_201805_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201805/skim_data_201805_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201805/skim_data_201805_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018May_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018May_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201805/skim_data_201805_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201805/skim_data_201805_part2_002_010.root");

//Jun 2018
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201806/skim_data_201806_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201806/skim_data_201806_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201806/skim_data_201806_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Jun_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Jun_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201806/skim_data_201806_part1_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201806/skim_data_201806_part2_002_010.root");

//Jul 2018
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201807/skim_data_201807_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201807/skim_data_201807_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201807/skim_data_201807_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Jul_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Jul_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201807/skim_data_201807_002_010_part1.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201807/skim_data_201807_002_010_part2.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201807/skim_data_201807_002_010_part3.root");

//Aug 2018
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201808/skim_data_201808_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201808/skim_data_201808_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201808/skim_data_201808_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Aug_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Aug_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201808/skim_data_201808_002_010_part1.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201808/skim_data_201808_002_010_part2.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201808/skim_data_201808_002_010_part3.root");

//Sep 2018
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201809/skim_data_201809_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201809/skim_data_201809_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201809/skim_data_201809_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Sep_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Sep_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201809/skim_data_201809_002_010_part1.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201809/skim_data_201809_002_010_part2.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201809/skim_data_201809_002_010_part3.root");

//Oct 2018
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201810/skim_data_201810_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201810/skim_data_201810_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201810/skim_data_201810_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Oct_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Oct_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201810/skim_data_201810_002_010_part1.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201810/skim_data_201810_002_010_part2.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201810/skim_data_201810_002_010_part3.root");

//Nov 2018
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201811/skim_data_201811_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201811/skim_data_201811_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201811/skim_data_201811_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Nov_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Nov_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201811/skim_data_201811_002_010_part1.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201811/skim_data_201811_002_010_part2.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201811/skim_data_201811_002_010_part3.root");

//Dec 2018
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201812/skim_data_201812_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201812/skim_data_201812_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201812/skim_data_201812_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Dec_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2018Dec_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201812/skim_data_201812_002_010_part1.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201812/skim_data_201812_002_010_part2.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201812/skim_data_201812_002_010_part3.root");

//Jan 2019
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201901/skim_data_201901_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201901/skim_data_201901_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201901/skim_data_201901_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Jan_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Jan_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201901/skim_data_201901_002_010_part1.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201901/skim_data_201901_002_010_part2.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201901/skim_data_201901_002_010_part3.root");

//Feb 2019
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201902/skim_data_201902_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201902/skim_data_201902_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201902/skim_data_201902_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Feb_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Feb_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201902/skim_data_201902_002_010_part1.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201902/skim_data_201902_002_010_part2.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201902/skim_data_201902_002_010_part3.root");

//Mar 2019
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201903/skim_data_201903_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201903/skim_data_201903_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201903/skim_data_201903_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Mar_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Mar_500_000_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201903/skim_data_201903_002_010_part1.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201903/skim_data_201903_002_010_part2.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201903/skim_data_201903_002_010_part3.root");

//Apr 2019
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201904/skim_data_201904_002_010_part1.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201904/skim_data_201904_002_010_part2.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201904/skim_data_201904_002_010_part3.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201904/skim_data_201904_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201904/skim_data_201904_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201904/skim_data_201904_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Apr_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Apr_500_000_SATcor_Protons.root");

//May 2019
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201905/skim_data_201905_002_010_part1.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201905/skim_data_201905_002_010_part2.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201905/skim_data_201905_002_010_part3.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201905/skim_data_201905_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201905/skim_data_201905_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201905/skim_data_201905_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019May_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019May_500_000_SATcor_Protons.root");

//Jun 2019
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201906/skim_data_201906_002_010_part1.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201906/skim_data_201906_002_010_part2.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201906/skim_data_201906_002_010_part3.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201906/skim_data_201906_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201906/skim_data_201906_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201906/skim_data_201906_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Jun_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Jun_500_000_SATcor_Protons.root");

//Jul 2019
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201907/skim_data_201907_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201907/skim_data_201907_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201907/skim_data_201907_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201907/skim_data_201907_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Jul_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Jul_500_000_SATcor_Protons.root");

//Ago 2019
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201908/skim_data_201908_002_010.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201908/skim_data_201908_010_025.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201908/skim_data_201908_025_050.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW/201908/skim_data_201908_050_100.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Aug_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Aug_500_000_SATcor_Protons.root");

//Sep 2019
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201909/skim_data_2019Sep_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201909/skim_data_2019Sep_010_025_SATcor_new.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201909/skim_data_2019Sep_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201909/skim_data_2019Sep_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Sep_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Sep_500_000_SATcor_Protons.root");

//Oct 2019
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201910/skim_data_2019Oct_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201910/skim_data_2019Oct_010_025_SATcor_new.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201910/skim_data_2019Oct_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201910/skim_data_2019Oct_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Oct_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Oct_500_000_SATcor_Protons.root");

//Nov 2019
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201911/skim_data_2019Nov_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201911/skim_data_2019Nov_010_025_SATcor_new.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201911/skim_data_2019Nov_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201911/skim_data_2019Nov_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Nov_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Nov_500_000_SATcor_Protons.root");

//Dec 2019
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201912/skim_data_2019Dec_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201912/skim_data_2019Dec_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201912/skim_data_2019Dec_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/201912/skim_data_2019Dec_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Dec_100_500_SATcor_Protons.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2019Dec_500_000_SATcor_Protons.root");

//Jan 2020
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202001/skim_data_2020Jan_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202001/skim_data_2020Jan_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202001/skim_data_2020Jan_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202001/skim_data_2020Jan_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Jan_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Jan_500_000_SATcor.root");

//Feb 2020
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202002/skim_data_2020Feb_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202002/skim_data_2020Feb_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202002/skim_data_2020Feb_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202002/skim_data_2020Feb_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Feb_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Feb_500_000_SATcor.root");

//Mar 2020
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202003/skim_data_2020Mar_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202003/skim_data_2020Mar_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202003/skim_data_2020Mar_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202003/skim_data_2020Mar_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Mar_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Mar_500_000_SATcor.root");

//Apr 2020
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202004/skim_data_2020Apr_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202004/skim_data_2020Apr_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202004/skim_data_2020Apr_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202004/skim_data_2020Apr_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Apr_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Apr_500_000_SATcor.root");

//May 2020
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202005/skim_data_2020May_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202005/skim_data_2020May_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202005/skim_data_2020May_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202005/skim_data_2020May_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020May_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020May_500_000_SATcor.root");

//Jun 2020
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202006/skim_data_2020Jun_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202006/skim_data_2020Jun_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202006/skim_data_2020Jun_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202006/skim_data_2020Jun_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Jun_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Jun_500_000_SATcor.root");

//Jul 2020
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202007/skim_data_2020Jul_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202007/skim_data_2020Jul_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202007/skim_data_2020Jul_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202007/skim_data_2020Jul_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Jul_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Jul_500_000_SATcor.root");

//Aug 2020
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202008/skim_data_2020Aug_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202008/skim_data_2020Aug_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202008/skim_data_2020Aug_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202008/skim_data_2020Aug_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Aug_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Aug_500_000_SATcor.root");

//Sep 2020
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202009/skim_data_2020Sep_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202009/skim_data_2020Sep_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202009/skim_data_2020Sep_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202009/skim_data_2020Sep_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Sep_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Sep_500_000_SATcor.root");

//Oct 2020
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202010/skim_data_2020Oct_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202010/skim_data_2020Oct_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202010/skim_data_2020Oct_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202010/skim_data_2020Oct_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Oct_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Oct_500_000_SATcor.root");

//Nov 2020
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202011/skim_data_2020Nov_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202011/skim_data_2020Nov_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202011/skim_data_2020Nov_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202011/skim_data_2020Nov_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Nov_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Nov_500_000_SATcor.root");

//Dec 2020
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202012/skim_data_2020Dec_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202012/skim_data_2020Dec_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202012/skim_data_2020Dec_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/FlightData/Skim_data_v600_v4_NEW_SATcor/202012/skim_data_2020Dec_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Dec_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2020Dec_500_000_SATcor.root");

//Jan 2021
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jan_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jan_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jan_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jan_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jan_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jan_500_000_SATcor.root");

//Feb 2021
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Feb_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Feb_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Feb_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Feb_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Feb_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Feb_500_000_SATcor.root");

//Mar 2021
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Mar_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Mar_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Mar_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Mar_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Mar_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Mar_500_000_SATcor.root");

//Apr 2021
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Apr_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Apr_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Apr_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Apr_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Apr_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Apr_500_000_SATcor.root");

//May 2021
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021May_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021May_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021May_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021May_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021May_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021May_500_000_SATcor.root");

//Jun 2021
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jun_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jun_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jun_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jun_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jun_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jun_500_000_SATcor.root");

//Jul 2021
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jul_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jul_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jul_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jul_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jul_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Jul_500_000_SATcor.root");

//Aug 2021
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Aug_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Aug_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Aug_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Aug_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Aug_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Aug_500_000_SATcor.root");

//Sep 2021
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Sep_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Sep_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Sep_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Sep_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Sep_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Sep_500_000_SATcor.root");

//Oct 2021
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Oct_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Oct_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Oct_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Oct_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Oct_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Oct_500_000_SATcor.root");

//Nov 2021
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Nov_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Nov_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Nov_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Nov_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Nov_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Nov_500_000_SATcor.root");

//Dec 2021
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Dec_002_010_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Dec_010_025_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Dec_025_050_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Dec_050_100_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Dec_100_500_SATcor.root");
skim->Add("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/FLIGHT/skim_data_2021Dec_500_000_SATcor.root");

cout<<"Orbital Data Enties:"<<skim -> GetEntries()<<endl;




TCut wcomP = "(2.*3.14159*3.14159*1.38*1.38)";
TCut wcomHe = "(2.*3.14159*3.14159)";
TCut GeoCorr = "(2.*3.14159*3.14159)/(2.*3.14159*3.14159*1.38*1.38)";
// 


TCut Trig_HEP="BGO_HET>0.";

TCut Trig_UnB="BGO_UBT>0.";
TCut Trig_AND="(BGO_HET>0.) && (BGO_UBT>0.)";
TCut Trig_UnBE="BGO_UBE>0.";

//
TCut cc204s = "(BGO_EnergyG>20.)";
TCut cc104s = "(BGO_EnergyG>10.)";

TCut cut00s = cc104s*Trig_HEP;

TCut cTrig="(BGO_EnergyG<1.e+7)"*Trig_HEP; 

//TCut Id_P = "(int(pvp_pdg)==2212)";


//
// --------------------------------------------------------------------------
// *** Selection CUTS ***
//
//
TCut cut00  = cc204s*Trig_HEP;
//
TCut cut01="PSD_psdY_Emax_layer>0.0 && PSD_psdX_Emax_layer>0.0";              // un segnale su entrambe i due piani di PSD
TCut cut02="BGO_EneLay[0]+BGO_EneLay[1]<BGO_EneLay[2]+BGO_EneLay[3]";

//TCut cut03="Etrack01>0.0 && Etrack23>0.0";                         // un segnale su entrambe i due piani di PSD
TCut cut03ups="(10.*(0.5*(TMath::Sign(1.,PSD_EnergyY_corr[0])+1)*PSD_EnergyY_corr[0]+0.5*(TMath::Sign(1.,PSD_EnergyY_corr[1])+1)*PSD_EnergyY_corr[1])/(0.5*(TMath::Sign(1.,PSD_EnergyY_corr[0])+1)*PSD_psdY_pathlength[0]+0.5*(TMath::Sign(1.,PSD_EnergyY_corr[1])+1)*PSD_psdY_pathlength[1]))>0."; // un segnale su entrambe i due piani di PSD
TCut cut03dws="(10.*(0.5*(TMath::Sign(1.,PSD_EnergyX_corr[0])+1)*PSD_EnergyX_corr[0]+0.5*(TMath::Sign(1.,PSD_EnergyX_corr[1])+1)*PSD_EnergyX_corr[1])/(0.5*(TMath::Sign(1.,PSD_EnergyX_corr[0])+1)*PSD_psdX_pathlength[0]+0.5*(TMath::Sign(1.,PSD_EnergyX_corr[1])+1)*PSD_psdX_pathlength[1]))>0."; // un segnale su entrambe i due piani di PSD
TCut cut03=cut03ups*cut03dws; 

//TCut cut04="(fabs(sqrt(Etrack01)-sqrt(epsdmax01/facslope))<5. && fabs(sqrt(Etrack23)-sqrt(epsdmax23/facslope))<5.)";
TCut cut04ups="abs(sqrt(10.*(0.5*(TMath::Sign(1.,PSD_EnergyY_corr[0])+1)*PSD_EnergyY_corr[0]+0.5*(TMath::Sign(1.,PSD_EnergyY_corr[1])+1)*PSD_EnergyY_corr[1])/(0.5*(TMath::Sign(1.,PSD_EnergyY_corr[0])+1)*PSD_psdY_pathlength[0]+0.5*(TMath::Sign(1.,PSD_EnergyY_corr[1])+1)*PSD_psdY_pathlength[1]))-sqrt(PSD_psdY_Emax_layer*STK_theta_correction))<5.";
TCut cut04dws="abs(sqrt(10.*(0.5*(TMath::Sign(1.,PSD_EnergyX_corr[0])+1)*PSD_EnergyX_corr[0]+0.5*(TMath::Sign(1.,PSD_EnergyX_corr[1])+1)*PSD_EnergyX_corr[1])/(0.5*(TMath::Sign(1.,PSD_EnergyX_corr[0])+1)*PSD_psdX_pathlength[0]+0.5*(TMath::Sign(1.,PSD_EnergyX_corr[1])+1)*PSD_psdX_pathlength[1]))-sqrt(PSD_psdX_Emax_layer*STK_theta_correction))<5.";
TCut cut04=cut04ups*cut04dws;


//TCut cut05="(fabs(cbgomax01-cbgostk01)<30.) && (fabs(cbgomax02-cbgostk02)<30.)";
TCut cut05="TMath::Abs(BGO_cbgomax[0]-BGO_cbgostk[0])<30. && TMath::Abs(BGO_cbgomax[1]-BGO_cbgostk[1])<30.";


 TCut ctot=cut00*cut01*cut02*cut03*cut04*cut05;

//
//TCut Contained="(ext==0)";
TCut Contained="1.";
// ............................ dal MC(p4)-Cont-LanGaus (-2,+3 sigma):
//TCut cpro1="(0.5*(Etrack01+Etrack23)-(0.008454*log10(Ebgotot*0.001)**4.036+2.111))<3*sqrt((0.007935*log10(Ebgotot*0.001)**3.432+0.1714)**2+0.085**2)";
//TCut cpro2="((0.008454*log10(Ebgotot*0.001)**4.036+2.111)-0.5*(Etrack01+Etrack23))<2*sqrt((0.007935*log10(Ebgotot*0.001)**3.432+0.1714)**2+0.085**2)";
//TCut cutP = cpro1*cpro2;


//************** NEW CHARGE SELECTION 06/04/2022 ***********************//

//TString HeMcMPVf="((15.78)+(-12.1)*log10(BGO_EnergyG)+(7.458)*log10(BGO_EnergyG)**2+(-1.978)*log10(BGO_EnergyG)**3+(0.2196)*log10(BGO_EnergyG)**4)";
//TString HeMcMPVf_Cor="((8.391)+(-0.03673)*log10(BGO_EnergyG)+(0.2714)*log10(BGO_EnergyG)**2+(-0.1093)*log10(BGO_EnergyG)**3+(0.02177)*log10(BGO_EnergyG)**4)";
TString HeMPVf="((7.852)+(1.035)*log10(BGO_EnergyG)+(-0.4836)*log10(BGO_EnergyG)**2+(0.1124)*log10(BGO_EnergyG)**3+(-0.001648)*log10(BGO_EnergyG)**4)";
//TString HeMcWidthf="((3.38)+(-4.905)*log10(BGO_EnergyG)+(3.076)*log10(BGO_EnergyG)**2+(-0.8565)*log10(BGO_EnergyG)**3+(0.09808)*log10(BGO_EnergyG)**4)";
//TString HeMcWidthf_Cor="((-0.1182)+(0.3145)*log10(BGO_EnergyG)+(0.01086)*log10(BGO_EnergyG)**2+(-0.02211)*log10(BGO_EnergyG)**3+(0.007157)*log10(BGO_EnergyG)**4)";
TString HeWidthf="((0.2505)+(0.02188)*log10(BGO_EnergyG)+(0.1554)*log10(BGO_EnergyG)**2+(-0.06573)*log10(BGO_EnergyG)**3+(0.01144)*log10(BGO_EnergyG)**4)";


//TString HeMcGSigmaf="(0.491)";
//TString HeMcGSigmaf_Cor="(0.02013)+(-4.928e-31)*log10(BGO_EnergyG)";
TString HeGSigmaf="(0.4789)+(-4.384e-33)*log10(BGO_EnergyG)";	


//TString HeMcFSig= "sqrt("+HeMcWidthf+"**2+"+HeMcGSigmaf+"**2)";
//TString HeMcFSigCor= "sqrt("+HeMcWidthf_Cor+"**2+"+HeMcGSigmaf_Cor+"**2)";
TString HeFSig= "sqrt("+HeWidthf+"**2+"+HeGSigmaf+"**2)";


//TString PMcMPVf="((4.576)+(-4.756)*log10(BGO_EnergyG)+(3.392)*log10(BGO_EnergyG)**2+(-1.036)*log10(BGO_EnergyG)**3+(0.129)*log10(BGO_EnergyG)**4)";
//TString PMcMPVf_Cor="((3.575)+(-3.177)*log10(BGO_EnergyG)+(2.398)*log10(BGO_EnergyG)**2+(-0.7764)*log10(BGO_EnergyG)**3+(0.09647)*log10(BGO_EnergyG)**4)";
TString PMPVf="((2.071)+(-0.3114)*log10(BGO_EnergyG)+(0.4153)*log10(BGO_EnergyG)**2+(-0.1876)*log10(BGO_EnergyG)**3+(0.03373)*log10(BGO_EnergyG)**4)";
//TString PMcWidthf="((1.112)+(-1.905)*log10(BGO_EnergyG)+(1.423)*log10(BGO_EnergyG)**2+(-0.4416)*log10(BGO_EnergyG)**3+(0.05443)*log10(BGO_EnergyG)**4)";
//TString PMcWidthf_Cor="((1.154)+(-2.034)*log10(BGO_EnergyG)+(1.513)*log10(BGO_EnergyG)**2+(-0.4764)*log10(BGO_EnergyG)**3+(0.05621)*log10(BGO_EnergyG)**4)";
TString PWidthf="((0.03606)+(0.09163)*log10(BGO_EnergyG)+(0.05757)*log10(BGO_EnergyG)**2+(-0.04983)*log10(BGO_EnergyG)**3+(0.01135)*log10(BGO_EnergyG)**4)";

//TString PMcGSigmaf="(1.0E-06)";
//TString PMcGSigmaf_Cor="(1.344E-09)+(-1.461E-40)*log10(BGO_EnergyG)";
TString PGSigmaf="(6.867E-09)+(9.713E-15)*log10(BGO_EnergyG)";		


//********************** END NEW **********************************//

//TString PMcFSig= "sqrt("+PMcWidthf+"**2+"+PMcGSigmaf+"**2)";
//TString PMcFSig_Cor= "sqrt("+PMcWidthf_Cor+"**2+"+PMcGSigmaf_Cor+"**2)";
TString PFSig= "sqrt("+PWidthf+"**2+"+PGSigmaf+"**2)";



TString p_low= "("+PMPVf+"-"+"(PSD_Etrack01+PSD_Etrack23)/2)<(2*"+PFSig+")";
TString he_high6= "((PSD_Etrack01+PSD_Etrack23)/2-"+HeMPVf+")<(6*"+HeFSig+")";

TCut PLow = p_low;
TCut HeHigh6= he_high6;


TCut Cut_PHe6 =PLow*HeHigh6;



TCut SpCut1= "log10(BGO_nhits_layers[12]*BGO_EneLay[12]/BGO_EnergyG+BGO_nhits_layers[13]*BGO_EneLay[13]/BGO_EnergyG)>-0.0001*(BGO_nbars)**2.0-1.34";
TCut SpCut2= "BGO_xtr > 16."; // nuovo taglio per gli elettroni

TCut SpCut=SpCut1*SpCut2;

TCut cut_00s  = cc104s*Trig_UnB;
TCut cut00_ANDs  = cc104s*Trig_HEP*Trig_UnB;

TCut ctotContUnB = cut_00s*cut01*cut02*cut03*cut04*cut05*SpCut*Contained*Cut_PHe6; 
TCut ctotContAND = cut00_ANDs*cut01*cut02*cut03*cut04*cut05*SpCut*Contained*Cut_PHe6; 
TCut ctotContUnBE = cc204s*Trig_UnBE*cut01*cut02*cut03*cut04*cut05*SpCut*Contained*Cut_PHe6; 
//TCut ctotCont6=Trig_UnBE*cut00*cut01*cut02*cut03*cut04*cut05*SpCut*Contained*Cut_PHe6;
//
TCut ctotCont=Trig_HEP*cut00*cut01*cut02*cut03*cut04*cut05*SpCut*Contained;
TCut ctotCont6=cut00*cut01*cut02*cut03*cut04*cut05*SpCut*Contained*Cut_PHe6;

// -----------------------------------------------------------------------------------------------------
//

//
// ----------------------------

TCut wspectP27 = "(MC_EnergyT)**(-1.7)";
TCut wspectHe27 = "(MC_EnergyT)**(-1.7)";



// for spectrum

int nbd = 48;
int ndec = 6;

 float arg1 = 1./float(nbd);
 int noe = nbd*ndec;
 cout << "Number of energy bin " << noe << endl;
 float e0 = 1.;
// float e0 = 2.;//shifted
 float Ebin[noe+1];
 Ebin[0] = e0;
 for (int i=1; i<noe+1; i++) {
    Ebin[i] = Ebin[i-1]*TMath::Power(10., arg1);
 }

cout<<"Now on Orbital Data..."<<endl;


// for spectrum
 
TH1F *h1SelBGO_orb6 = new TH1F("h1SelBGO_orb6", "Selected(E_bgo) orbital", noe, Ebin);//Orbital Data
/*
TH1F *h1Sel_UnB = new TH1F("h1Sel_UnB", "Selected(E_bgo) orbital", noe, Ebin);
TH1F *h1Sel_UnBE = new TH1F("h1Sel_UnBE", "Selected(E_bgo) orbital", noe, Ebin);
TH1F *h1Sel_AND = new TH1F("h1Sel_AND", "Selected(E_bgo) orbital", noe, Ebin);
*/
h1SelBGO_orb6->Sumw2();
/*
h1Sel_UnB->Sumw2();
h1Sel_UnBE->Sumw2();
h1Sel_AND->Sumw2();
*/
skim->Draw("(BGO_EnergyG)>>+h1SelBGO_orb6",ctotCont6,""); 
/*
skim->Draw("(BGO_EnergyG)>>+h1Sel_UnB",ctotContUnB,""); 
skim->Draw("(BGO_EnergyG)>>+h1Sel_UnBE",ctotContUnBE,""); 
skim->Draw("(BGO_EnergyG)>>+h1Sel_AND",ctotContAND,""); 
*/
cout<<"6sigma finished..."<<endl;

TFile *fout1 = new TFile("PHe_skim_Orb72Month_He6_48binperdecadeTEST.root","RECREATE");

// for spectrum

fout1->cd();

h1SelBGO_orb6->Write();
/*
h1Sel_UnB->Write();
h1Sel_UnBE->Write();
h1Sel_AND->Write();
*/
fout1->Close();

cout<<"End script."<<endl;
};

