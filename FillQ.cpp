// Automatically generated file using genQ.r script, don't modify!

#include "TransitionMatrix.h"

double TransitionMatrix::fillQ(double aOmega, double aK)
{
	mQ[   0] = -aK*mCodonFreq[1]-aOmega*mCodonFreq[2]-aOmega*mCodonFreq[3]-aK*aOmega*mCodonFreq[4]-aOmega*mCodonFreq[8]-aOmega*mCodonFreq[10]-aK*aOmega*mCodonFreq[13]-aOmega*mCodonFreq[29]-aOmega*mCodonFreq[45];
	mQ[   1] = aK*mCodonFreq[1];
	mQ[   2] = aOmega*mCodonFreq[2];
	mQ[   3] = aOmega*mCodonFreq[3];
	mQ[   4] = aK*aOmega*mCodonFreq[4];
	mQ[   8] = aOmega*mCodonFreq[8];
	mQ[  10] = aOmega*mCodonFreq[10];
	mQ[  13] = aK*aOmega*mCodonFreq[13];
	mQ[  29] = aOmega*mCodonFreq[29];
	mQ[  45] = aOmega*mCodonFreq[45];
	mQ[  61] = aK*mCodonFreq[0];
	mQ[  62] = -aK*mCodonFreq[0]-aOmega*mCodonFreq[2]-aOmega*mCodonFreq[3]-aK*aOmega*mCodonFreq[5]-aOmega*mCodonFreq[9]-aOmega*mCodonFreq[11]-aK*aOmega*mCodonFreq[14]-aOmega*mCodonFreq[30]-aOmega*mCodonFreq[46];
	mQ[  63] = aOmega*mCodonFreq[2];
	mQ[  64] = aOmega*mCodonFreq[3];
	mQ[  66] = aK*aOmega*mCodonFreq[5];
	mQ[  70] = aOmega*mCodonFreq[9];
	mQ[  72] = aOmega*mCodonFreq[11];
	mQ[  75] = aK*aOmega*mCodonFreq[14];
	mQ[  91] = aOmega*mCodonFreq[30];
	mQ[ 107] = aOmega*mCodonFreq[46];
	mQ[ 122] = aOmega*mCodonFreq[0];
	mQ[ 123] = aOmega*mCodonFreq[1];
	mQ[ 124] = -aOmega*mCodonFreq[0]-aOmega*mCodonFreq[1]-aK*mCodonFreq[3]-aK*aOmega*mCodonFreq[6]-aK*mCodonFreq[15]-aOmega*mCodonFreq[31]-aOmega*mCodonFreq[47];
	mQ[ 125] = aK*mCodonFreq[3];
	mQ[ 128] = aK*aOmega*mCodonFreq[6];
	mQ[ 137] = aK*mCodonFreq[15];
	mQ[ 153] = aOmega*mCodonFreq[31];
	mQ[ 169] = aOmega*mCodonFreq[47];
	mQ[ 183] = aOmega*mCodonFreq[0];
	mQ[ 184] = aOmega*mCodonFreq[1];
	mQ[ 185] = aK*mCodonFreq[2];
	mQ[ 186] = -aOmega*mCodonFreq[0]-aOmega*mCodonFreq[1]-aK*mCodonFreq[2]-aK*aOmega*mCodonFreq[7]-aOmega*mCodonFreq[12]-aK*mCodonFreq[16]-aOmega*mCodonFreq[32]-aOmega*mCodonFreq[48];
	mQ[ 190] = aK*aOmega*mCodonFreq[7];
	mQ[ 195] = aOmega*mCodonFreq[12];
	mQ[ 199] = aK*mCodonFreq[16];
	mQ[ 215] = aOmega*mCodonFreq[32];
	mQ[ 231] = aOmega*mCodonFreq[48];
	mQ[ 244] = aK*aOmega*mCodonFreq[0];
	mQ[ 248] = -aK*aOmega*mCodonFreq[0]-aK*mCodonFreq[5]-mCodonFreq[6]-mCodonFreq[7]-aOmega*mCodonFreq[8]-aOmega*mCodonFreq[10]-aK*aOmega*mCodonFreq[17]-aOmega*mCodonFreq[33]-aOmega*mCodonFreq[49];
	mQ[ 249] = aK*mCodonFreq[5];
	mQ[ 250] = mCodonFreq[6];
	mQ[ 251] = mCodonFreq[7];
	mQ[ 252] = aOmega*mCodonFreq[8];
	mQ[ 254] = aOmega*mCodonFreq[10];
	mQ[ 261] = aK*aOmega*mCodonFreq[17];
	mQ[ 277] = aOmega*mCodonFreq[33];
	mQ[ 293] = aOmega*mCodonFreq[49];
	mQ[ 306] = aK*aOmega*mCodonFreq[1];
	mQ[ 309] = aK*mCodonFreq[4];
	mQ[ 310] = -aK*aOmega*mCodonFreq[1]-aK*mCodonFreq[4]-mCodonFreq[6]-mCodonFreq[7]-aOmega*mCodonFreq[9]-aOmega*mCodonFreq[11]-aK*aOmega*mCodonFreq[18]-aOmega*mCodonFreq[34]-aOmega*mCodonFreq[50];
	mQ[ 311] = mCodonFreq[6];
	mQ[ 312] = mCodonFreq[7];
	mQ[ 314] = aOmega*mCodonFreq[9];
	mQ[ 316] = aOmega*mCodonFreq[11];
	mQ[ 323] = aK*aOmega*mCodonFreq[18];
	mQ[ 339] = aOmega*mCodonFreq[34];
	mQ[ 355] = aOmega*mCodonFreq[50];
	mQ[ 368] = aK*aOmega*mCodonFreq[2];
	mQ[ 370] = mCodonFreq[4];
	mQ[ 371] = mCodonFreq[5];
	mQ[ 372] = -aK*aOmega*mCodonFreq[2]-mCodonFreq[4]-mCodonFreq[5]-aK*mCodonFreq[7]-aK*aOmega*mCodonFreq[19]-aOmega*mCodonFreq[35]-aOmega*mCodonFreq[51];
	mQ[ 373] = aK*mCodonFreq[7];
	mQ[ 385] = aK*aOmega*mCodonFreq[19];
	mQ[ 401] = aOmega*mCodonFreq[35];
	mQ[ 417] = aOmega*mCodonFreq[51];
	mQ[ 430] = aK*aOmega*mCodonFreq[3];
	mQ[ 431] = mCodonFreq[4];
	mQ[ 432] = mCodonFreq[5];
	mQ[ 433] = aK*mCodonFreq[6];
	mQ[ 434] = -aK*aOmega*mCodonFreq[3]-mCodonFreq[4]-mCodonFreq[5]-aK*mCodonFreq[6]-aOmega*mCodonFreq[12]-aK*aOmega*mCodonFreq[20]-aOmega*mCodonFreq[36]-aOmega*mCodonFreq[52];
	mQ[ 439] = aOmega*mCodonFreq[12];
	mQ[ 447] = aK*aOmega*mCodonFreq[20];
	mQ[ 463] = aOmega*mCodonFreq[36];
	mQ[ 479] = aOmega*mCodonFreq[52];
	mQ[ 488] = aOmega*mCodonFreq[0];
	mQ[ 492] = aOmega*mCodonFreq[4];
	mQ[ 496] = -aOmega*mCodonFreq[0]-aOmega*mCodonFreq[4]-aK*mCodonFreq[9]-aK*aOmega*mCodonFreq[10]-aK*aOmega*mCodonFreq[21]-aOmega*mCodonFreq[37]-aOmega*mCodonFreq[53];
	mQ[ 497] = aK*mCodonFreq[9];
	mQ[ 498] = aK*aOmega*mCodonFreq[10];
	mQ[ 509] = aK*aOmega*mCodonFreq[21];
	mQ[ 525] = aOmega*mCodonFreq[37];
	mQ[ 541] = aOmega*mCodonFreq[53];
	mQ[ 550] = aOmega*mCodonFreq[1];
	mQ[ 554] = aOmega*mCodonFreq[5];
	mQ[ 557] = aK*mCodonFreq[8];
	mQ[ 558] = -aOmega*mCodonFreq[1]-aOmega*mCodonFreq[5]-aK*mCodonFreq[8]-aK*aOmega*mCodonFreq[11]-aK*aOmega*mCodonFreq[22]-aOmega*mCodonFreq[38]-aOmega*mCodonFreq[54];
	mQ[ 560] = aK*aOmega*mCodonFreq[11];
	mQ[ 571] = aK*aOmega*mCodonFreq[22];
	mQ[ 587] = aOmega*mCodonFreq[38];
	mQ[ 603] = aOmega*mCodonFreq[54];
	mQ[ 610] = aOmega*mCodonFreq[0];
	mQ[ 614] = aOmega*mCodonFreq[4];
	mQ[ 618] = aK*aOmega*mCodonFreq[8];
	mQ[ 620] = -aOmega*mCodonFreq[0]-aOmega*mCodonFreq[4]-aK*aOmega*mCodonFreq[8]-aK*mCodonFreq[11]-aOmega*mCodonFreq[12]-aK*aOmega*mCodonFreq[25]-aOmega*mCodonFreq[41]-aOmega*mCodonFreq[57];
	mQ[ 621] = aK*mCodonFreq[11];
	mQ[ 622] = aOmega*mCodonFreq[12];
	mQ[ 635] = aK*aOmega*mCodonFreq[25];
	mQ[ 651] = aOmega*mCodonFreq[41];
	mQ[ 667] = aOmega*mCodonFreq[57];
	mQ[ 672] = aOmega*mCodonFreq[1];
	mQ[ 676] = aOmega*mCodonFreq[5];
	mQ[ 680] = aK*aOmega*mCodonFreq[9];
	mQ[ 681] = aK*mCodonFreq[10];
	mQ[ 682] = -aOmega*mCodonFreq[1]-aOmega*mCodonFreq[5]-aK*aOmega*mCodonFreq[9]-aK*mCodonFreq[10]-aOmega*mCodonFreq[12]-aK*aOmega*mCodonFreq[26]-aOmega*mCodonFreq[42]-aOmega*mCodonFreq[58];
	mQ[ 683] = aOmega*mCodonFreq[12];
	mQ[ 697] = aK*aOmega*mCodonFreq[26];
	mQ[ 713] = aOmega*mCodonFreq[42];
	mQ[ 729] = aOmega*mCodonFreq[58];
	mQ[ 735] = aOmega*mCodonFreq[3];
	mQ[ 739] = aOmega*mCodonFreq[7];
	mQ[ 742] = aOmega*mCodonFreq[10];
	mQ[ 743] = aOmega*mCodonFreq[11];
	mQ[ 744] = -aOmega*mCodonFreq[3]-aOmega*mCodonFreq[7]-aOmega*mCodonFreq[10]-aOmega*mCodonFreq[11]-aK*aOmega*mCodonFreq[28]-aOmega*mCodonFreq[44]-aOmega*mCodonFreq[60];
	mQ[ 760] = aK*aOmega*mCodonFreq[28];
	mQ[ 776] = aOmega*mCodonFreq[44];
	mQ[ 792] = aOmega*mCodonFreq[60];
	mQ[ 793] = aK*aOmega*mCodonFreq[0];
	mQ[ 806] = -aK*aOmega*mCodonFreq[0]-aK*mCodonFreq[14]-mCodonFreq[15]-mCodonFreq[16]-aK*aOmega*mCodonFreq[17]-aOmega*mCodonFreq[21]-aOmega*mCodonFreq[25]-aOmega*mCodonFreq[29]-aOmega*mCodonFreq[45];
	mQ[ 807] = aK*mCodonFreq[14];
	mQ[ 808] = mCodonFreq[15];
	mQ[ 809] = mCodonFreq[16];
	mQ[ 810] = aK*aOmega*mCodonFreq[17];
	mQ[ 814] = aOmega*mCodonFreq[21];
	mQ[ 818] = aOmega*mCodonFreq[25];
	mQ[ 822] = aOmega*mCodonFreq[29];
	mQ[ 838] = aOmega*mCodonFreq[45];
	mQ[ 855] = aK*aOmega*mCodonFreq[1];
	mQ[ 867] = aK*mCodonFreq[13];
	mQ[ 868] = -aK*aOmega*mCodonFreq[1]-aK*mCodonFreq[13]-mCodonFreq[15]-mCodonFreq[16]-aK*aOmega*mCodonFreq[18]-aOmega*mCodonFreq[22]-aOmega*mCodonFreq[26]-aOmega*mCodonFreq[30]-aOmega*mCodonFreq[46];
	mQ[ 869] = mCodonFreq[15];
	mQ[ 870] = mCodonFreq[16];
	mQ[ 872] = aK*aOmega*mCodonFreq[18];
	mQ[ 876] = aOmega*mCodonFreq[22];
	mQ[ 880] = aOmega*mCodonFreq[26];
	mQ[ 884] = aOmega*mCodonFreq[30];
	mQ[ 900] = aOmega*mCodonFreq[46];
	mQ[ 917] = aK*mCodonFreq[2];
	mQ[ 928] = mCodonFreq[13];
	mQ[ 929] = mCodonFreq[14];
	mQ[ 930] = -aK*mCodonFreq[2]-mCodonFreq[13]-mCodonFreq[14]-aK*mCodonFreq[16]-aK*aOmega*mCodonFreq[19]-aOmega*mCodonFreq[23]-aOmega*mCodonFreq[27]-aOmega*mCodonFreq[31]-aOmega*mCodonFreq[47];
	mQ[ 931] = aK*mCodonFreq[16];
	mQ[ 934] = aK*aOmega*mCodonFreq[19];
	mQ[ 938] = aOmega*mCodonFreq[23];
	mQ[ 942] = aOmega*mCodonFreq[27];
	mQ[ 946] = aOmega*mCodonFreq[31];
	mQ[ 962] = aOmega*mCodonFreq[47];
	mQ[ 979] = aK*mCodonFreq[3];
	mQ[ 989] = mCodonFreq[13];
	mQ[ 990] = mCodonFreq[14];
	mQ[ 991] = aK*mCodonFreq[15];
	mQ[ 992] = -aK*mCodonFreq[3]-mCodonFreq[13]-mCodonFreq[14]-aK*mCodonFreq[15]-aK*aOmega*mCodonFreq[20]-aOmega*mCodonFreq[24]-aOmega*mCodonFreq[28]-aOmega*mCodonFreq[32]-aOmega*mCodonFreq[48];
	mQ[ 996] = aK*aOmega*mCodonFreq[20];
	mQ[1000] = aOmega*mCodonFreq[24];
	mQ[1004] = aOmega*mCodonFreq[28];
	mQ[1008] = aOmega*mCodonFreq[32];
	mQ[1024] = aOmega*mCodonFreq[48];
	mQ[1041] = aK*aOmega*mCodonFreq[4];
	mQ[1050] = aK*aOmega*mCodonFreq[13];
	mQ[1054] = -aK*aOmega*mCodonFreq[4]-aK*aOmega*mCodonFreq[13]-aK*mCodonFreq[18]-mCodonFreq[19]-mCodonFreq[20]-aOmega*mCodonFreq[21]-aOmega*mCodonFreq[25]-aOmega*mCodonFreq[33]-aOmega*mCodonFreq[49];
	mQ[1055] = aK*mCodonFreq[18];
	mQ[1056] = mCodonFreq[19];
	mQ[1057] = mCodonFreq[20];
	mQ[1058] = aOmega*mCodonFreq[21];
	mQ[1062] = aOmega*mCodonFreq[25];
	mQ[1070] = aOmega*mCodonFreq[33];
	mQ[1086] = aOmega*mCodonFreq[49];
	mQ[1103] = aK*aOmega*mCodonFreq[5];
	mQ[1112] = aK*aOmega*mCodonFreq[14];
	mQ[1115] = aK*mCodonFreq[17];
	mQ[1116] = -aK*aOmega*mCodonFreq[5]-aK*aOmega*mCodonFreq[14]-aK*mCodonFreq[17]-mCodonFreq[19]-mCodonFreq[20]-aOmega*mCodonFreq[22]-aOmega*mCodonFreq[26]-aOmega*mCodonFreq[34]-aOmega*mCodonFreq[50];
	mQ[1117] = mCodonFreq[19];
	mQ[1118] = mCodonFreq[20];
	mQ[1120] = aOmega*mCodonFreq[22];
	mQ[1124] = aOmega*mCodonFreq[26];
	mQ[1132] = aOmega*mCodonFreq[34];
	mQ[1148] = aOmega*mCodonFreq[50];
	mQ[1165] = aK*aOmega*mCodonFreq[6];
	mQ[1174] = aK*aOmega*mCodonFreq[15];
	mQ[1176] = mCodonFreq[17];
	mQ[1177] = mCodonFreq[18];
	mQ[1178] = -aK*aOmega*mCodonFreq[6]-aK*aOmega*mCodonFreq[15]-mCodonFreq[17]-mCodonFreq[18]-aK*mCodonFreq[20]-aOmega*mCodonFreq[23]-aOmega*mCodonFreq[27]-aOmega*mCodonFreq[35]-aOmega*mCodonFreq[51];
	mQ[1179] = aK*mCodonFreq[20];
	mQ[1182] = aOmega*mCodonFreq[23];
	mQ[1186] = aOmega*mCodonFreq[27];
	mQ[1194] = aOmega*mCodonFreq[35];
	mQ[1210] = aOmega*mCodonFreq[51];
	mQ[1227] = aK*aOmega*mCodonFreq[7];
	mQ[1236] = aK*aOmega*mCodonFreq[16];
	mQ[1237] = mCodonFreq[17];
	mQ[1238] = mCodonFreq[18];
	mQ[1239] = aK*mCodonFreq[19];
	mQ[1240] = -aK*aOmega*mCodonFreq[7]-aK*aOmega*mCodonFreq[16]-mCodonFreq[17]-mCodonFreq[18]-aK*mCodonFreq[19]-aOmega*mCodonFreq[24]-aOmega*mCodonFreq[28]-aOmega*mCodonFreq[36]-aOmega*mCodonFreq[52];
	mQ[1244] = aOmega*mCodonFreq[24];
	mQ[1248] = aOmega*mCodonFreq[28];
	mQ[1256] = aOmega*mCodonFreq[36];
	mQ[1272] = aOmega*mCodonFreq[52];
	mQ[1289] = aK*aOmega*mCodonFreq[8];
	mQ[1294] = aOmega*mCodonFreq[13];
	mQ[1298] = aOmega*mCodonFreq[17];
	mQ[1302] = -aK*aOmega*mCodonFreq[8]-aOmega*mCodonFreq[13]-aOmega*mCodonFreq[17]-aK*mCodonFreq[22]-aOmega*mCodonFreq[23]-aOmega*mCodonFreq[24]-aK*aOmega*mCodonFreq[25]-aOmega*mCodonFreq[37]-aOmega*mCodonFreq[53];
	mQ[1303] = aK*mCodonFreq[22];
	mQ[1304] = aOmega*mCodonFreq[23];
	mQ[1305] = aOmega*mCodonFreq[24];
	mQ[1306] = aK*aOmega*mCodonFreq[25];
	mQ[1318] = aOmega*mCodonFreq[37];
	mQ[1334] = aOmega*mCodonFreq[53];
	mQ[1351] = aK*aOmega*mCodonFreq[9];
	mQ[1356] = aOmega*mCodonFreq[14];
	mQ[1360] = aOmega*mCodonFreq[18];
	mQ[1363] = aK*mCodonFreq[21];
	mQ[1364] = -aK*aOmega*mCodonFreq[9]-aOmega*mCodonFreq[14]-aOmega*mCodonFreq[18]-aK*mCodonFreq[21]-aOmega*mCodonFreq[23]-aOmega*mCodonFreq[24]-aK*aOmega*mCodonFreq[26]-aOmega*mCodonFreq[38]-aOmega*mCodonFreq[54];
	mQ[1365] = aOmega*mCodonFreq[23];
	mQ[1366] = aOmega*mCodonFreq[24];
	mQ[1368] = aK*aOmega*mCodonFreq[26];
	mQ[1380] = aOmega*mCodonFreq[38];
	mQ[1396] = aOmega*mCodonFreq[54];
	mQ[1418] = aOmega*mCodonFreq[15];
	mQ[1422] = aOmega*mCodonFreq[19];
	mQ[1424] = aOmega*mCodonFreq[21];
	mQ[1425] = aOmega*mCodonFreq[22];
	mQ[1426] = -aOmega*mCodonFreq[15]-aOmega*mCodonFreq[19]-aOmega*mCodonFreq[21]-aOmega*mCodonFreq[22]-aK*mCodonFreq[24]-aK*aOmega*mCodonFreq[27]-aOmega*mCodonFreq[39]-aOmega*mCodonFreq[55];
	mQ[1427] = aK*mCodonFreq[24];
	mQ[1430] = aK*aOmega*mCodonFreq[27];
	mQ[1442] = aOmega*mCodonFreq[39];
	mQ[1458] = aOmega*mCodonFreq[55];
	mQ[1480] = aOmega*mCodonFreq[16];
	mQ[1484] = aOmega*mCodonFreq[20];
	mQ[1485] = aOmega*mCodonFreq[21];
	mQ[1486] = aOmega*mCodonFreq[22];
	mQ[1487] = aK*mCodonFreq[23];
	mQ[1488] = -aOmega*mCodonFreq[16]-aOmega*mCodonFreq[20]-aOmega*mCodonFreq[21]-aOmega*mCodonFreq[22]-aK*mCodonFreq[23]-aK*aOmega*mCodonFreq[28]-aOmega*mCodonFreq[40]-aOmega*mCodonFreq[56];
	mQ[1492] = aK*aOmega*mCodonFreq[28];
	mQ[1504] = aOmega*mCodonFreq[40];
	mQ[1520] = aOmega*mCodonFreq[56];
	mQ[1535] = aK*aOmega*mCodonFreq[10];
	mQ[1538] = aOmega*mCodonFreq[13];
	mQ[1542] = aOmega*mCodonFreq[17];
	mQ[1546] = aK*aOmega*mCodonFreq[21];
	mQ[1550] = -aK*aOmega*mCodonFreq[10]-aOmega*mCodonFreq[13]-aOmega*mCodonFreq[17]-aK*aOmega*mCodonFreq[21]-aK*mCodonFreq[26]-mCodonFreq[27]-mCodonFreq[28]-aOmega*mCodonFreq[41]-aOmega*mCodonFreq[57];
	mQ[1551] = aK*mCodonFreq[26];
	mQ[1552] = mCodonFreq[27];
	mQ[1553] = mCodonFreq[28];
	mQ[1566] = aOmega*mCodonFreq[41];
	mQ[1582] = aOmega*mCodonFreq[57];
	mQ[1597] = aK*aOmega*mCodonFreq[11];
	mQ[1600] = aOmega*mCodonFreq[14];
	mQ[1604] = aOmega*mCodonFreq[18];
	mQ[1608] = aK*aOmega*mCodonFreq[22];
	mQ[1611] = aK*mCodonFreq[25];
	mQ[1612] = -aK*aOmega*mCodonFreq[11]-aOmega*mCodonFreq[14]-aOmega*mCodonFreq[18]-aK*aOmega*mCodonFreq[22]-aK*mCodonFreq[25]-mCodonFreq[27]-mCodonFreq[28]-aOmega*mCodonFreq[42]-aOmega*mCodonFreq[58];
	mQ[1613] = mCodonFreq[27];
	mQ[1614] = mCodonFreq[28];
	mQ[1628] = aOmega*mCodonFreq[42];
	mQ[1644] = aOmega*mCodonFreq[58];
	mQ[1662] = aOmega*mCodonFreq[15];
	mQ[1666] = aOmega*mCodonFreq[19];
	mQ[1670] = aK*aOmega*mCodonFreq[23];
	mQ[1672] = mCodonFreq[25];
	mQ[1673] = mCodonFreq[26];
	mQ[1674] = -aOmega*mCodonFreq[15]-aOmega*mCodonFreq[19]-aK*aOmega*mCodonFreq[23]-mCodonFreq[25]-mCodonFreq[26]-aK*mCodonFreq[28]-mCodonFreq[43]-aOmega*mCodonFreq[59];
	mQ[1675] = aK*mCodonFreq[28];
	mQ[1690] = mCodonFreq[43];
	mQ[1706] = aOmega*mCodonFreq[59];
	mQ[1720] = aK*aOmega*mCodonFreq[12];
	mQ[1724] = aOmega*mCodonFreq[16];
	mQ[1728] = aOmega*mCodonFreq[20];
	mQ[1732] = aK*aOmega*mCodonFreq[24];
	mQ[1733] = mCodonFreq[25];
	mQ[1734] = mCodonFreq[26];
	mQ[1735] = aK*mCodonFreq[27];
	mQ[1736] = -aK*aOmega*mCodonFreq[12]-aOmega*mCodonFreq[16]-aOmega*mCodonFreq[20]-aK*aOmega*mCodonFreq[24]-mCodonFreq[25]-mCodonFreq[26]-aK*mCodonFreq[27]-mCodonFreq[44]-aOmega*mCodonFreq[60];
	mQ[1752] = mCodonFreq[44];
	mQ[1768] = aOmega*mCodonFreq[60];
	mQ[1769] = aOmega*mCodonFreq[0];
	mQ[1782] = aOmega*mCodonFreq[13];
	mQ[1798] = -aOmega*mCodonFreq[0]-aOmega*mCodonFreq[13]-aK*mCodonFreq[30]-mCodonFreq[31]-aOmega*mCodonFreq[32]-aK*aOmega*mCodonFreq[33]-aOmega*mCodonFreq[37]-aOmega*mCodonFreq[41]-aK*aOmega*mCodonFreq[45];
	mQ[1799] = aK*mCodonFreq[30];
	mQ[1800] = mCodonFreq[31];
	mQ[1801] = aOmega*mCodonFreq[32];
	mQ[1802] = aK*aOmega*mCodonFreq[33];
	mQ[1806] = aOmega*mCodonFreq[37];
	mQ[1810] = aOmega*mCodonFreq[41];
	mQ[1814] = aK*aOmega*mCodonFreq[45];
	mQ[1831] = aOmega*mCodonFreq[1];
	mQ[1844] = aOmega*mCodonFreq[14];
	mQ[1859] = aK*mCodonFreq[29];
	mQ[1860] = -aOmega*mCodonFreq[1]-aOmega*mCodonFreq[14]-aK*mCodonFreq[29]-mCodonFreq[31]-aOmega*mCodonFreq[32]-aK*aOmega*mCodonFreq[34]-aOmega*mCodonFreq[38]-aOmega*mCodonFreq[42]-aK*aOmega*mCodonFreq[46];
	mQ[1861] = mCodonFreq[31];
	mQ[1862] = aOmega*mCodonFreq[32];
	mQ[1864] = aK*aOmega*mCodonFreq[34];
	mQ[1868] = aOmega*mCodonFreq[38];
	mQ[1872] = aOmega*mCodonFreq[42];
	mQ[1876] = aK*aOmega*mCodonFreq[46];
	mQ[1893] = aOmega*mCodonFreq[2];
	mQ[1906] = aOmega*mCodonFreq[15];
	mQ[1920] = mCodonFreq[29];
	mQ[1921] = mCodonFreq[30];
	mQ[1922] = -aOmega*mCodonFreq[2]-aOmega*mCodonFreq[15]-mCodonFreq[29]-mCodonFreq[30]-aK*aOmega*mCodonFreq[32]-aK*aOmega*mCodonFreq[35]-aOmega*mCodonFreq[39]-aOmega*mCodonFreq[43]-aK*aOmega*mCodonFreq[47];
	mQ[1923] = aK*aOmega*mCodonFreq[32];
	mQ[1926] = aK*aOmega*mCodonFreq[35];
	mQ[1930] = aOmega*mCodonFreq[39];
	mQ[1934] = aOmega*mCodonFreq[43];
	mQ[1938] = aK*aOmega*mCodonFreq[47];
	mQ[1955] = aOmega*mCodonFreq[3];
	mQ[1968] = aOmega*mCodonFreq[16];
	mQ[1981] = aOmega*mCodonFreq[29];
	mQ[1982] = aOmega*mCodonFreq[30];
	mQ[1983] = aK*aOmega*mCodonFreq[31];
	mQ[1984] = -aOmega*mCodonFreq[3]-aOmega*mCodonFreq[16]-aOmega*mCodonFreq[29]-aOmega*mCodonFreq[30]-aK*aOmega*mCodonFreq[31]-aK*aOmega*mCodonFreq[36]-aOmega*mCodonFreq[40]-aOmega*mCodonFreq[44]-aK*aOmega*mCodonFreq[48];
	mQ[1988] = aK*aOmega*mCodonFreq[36];
	mQ[1992] = aOmega*mCodonFreq[40];
	mQ[1996] = aOmega*mCodonFreq[44];
	mQ[2000] = aK*aOmega*mCodonFreq[48];
	mQ[2017] = aOmega*mCodonFreq[4];
	mQ[2030] = aOmega*mCodonFreq[17];
	mQ[2042] = aK*aOmega*mCodonFreq[29];
	mQ[2046] = -aOmega*mCodonFreq[4]-aOmega*mCodonFreq[17]-aK*aOmega*mCodonFreq[29]-aK*mCodonFreq[34]-mCodonFreq[35]-mCodonFreq[36]-aOmega*mCodonFreq[37]-aOmega*mCodonFreq[41]-aK*aOmega*mCodonFreq[49];
	mQ[2047] = aK*mCodonFreq[34];
	mQ[2048] = mCodonFreq[35];
	mQ[2049] = mCodonFreq[36];
	mQ[2050] = aOmega*mCodonFreq[37];
	mQ[2054] = aOmega*mCodonFreq[41];
	mQ[2062] = aK*aOmega*mCodonFreq[49];
	mQ[2079] = aOmega*mCodonFreq[5];
	mQ[2092] = aOmega*mCodonFreq[18];
	mQ[2104] = aK*aOmega*mCodonFreq[30];
	mQ[2107] = aK*mCodonFreq[33];
	mQ[2108] = -aOmega*mCodonFreq[5]-aOmega*mCodonFreq[18]-aK*aOmega*mCodonFreq[30]-aK*mCodonFreq[33]-mCodonFreq[35]-mCodonFreq[36]-aOmega*mCodonFreq[38]-aOmega*mCodonFreq[42]-aK*aOmega*mCodonFreq[50];
	mQ[2109] = mCodonFreq[35];
	mQ[2110] = mCodonFreq[36];
	mQ[2112] = aOmega*mCodonFreq[38];
	mQ[2116] = aOmega*mCodonFreq[42];
	mQ[2124] = aK*aOmega*mCodonFreq[50];
	mQ[2141] = aOmega*mCodonFreq[6];
	mQ[2154] = aOmega*mCodonFreq[19];
	mQ[2166] = aK*aOmega*mCodonFreq[31];
	mQ[2168] = mCodonFreq[33];
	mQ[2169] = mCodonFreq[34];
	mQ[2170] = -aOmega*mCodonFreq[6]-aOmega*mCodonFreq[19]-aK*aOmega*mCodonFreq[31]-mCodonFreq[33]-mCodonFreq[34]-aK*mCodonFreq[36]-aOmega*mCodonFreq[39]-aOmega*mCodonFreq[43]-aK*aOmega*mCodonFreq[51];
	mQ[2171] = aK*mCodonFreq[36];
	mQ[2174] = aOmega*mCodonFreq[39];
	mQ[2178] = aOmega*mCodonFreq[43];
	mQ[2186] = aK*aOmega*mCodonFreq[51];
	mQ[2203] = aOmega*mCodonFreq[7];
	mQ[2216] = aOmega*mCodonFreq[20];
	mQ[2228] = aK*aOmega*mCodonFreq[32];
	mQ[2229] = mCodonFreq[33];
	mQ[2230] = mCodonFreq[34];
	mQ[2231] = aK*mCodonFreq[35];
	mQ[2232] = -aOmega*mCodonFreq[7]-aOmega*mCodonFreq[20]-aK*aOmega*mCodonFreq[32]-mCodonFreq[33]-mCodonFreq[34]-aK*mCodonFreq[35]-aOmega*mCodonFreq[40]-aOmega*mCodonFreq[44]-aK*aOmega*mCodonFreq[52];
	mQ[2236] = aOmega*mCodonFreq[40];
	mQ[2240] = aOmega*mCodonFreq[44];
	mQ[2248] = aK*aOmega*mCodonFreq[52];
	mQ[2265] = aOmega*mCodonFreq[8];
	mQ[2278] = aOmega*mCodonFreq[21];
	mQ[2286] = aOmega*mCodonFreq[29];
	mQ[2290] = aOmega*mCodonFreq[33];
	mQ[2294] = -aOmega*mCodonFreq[8]-aOmega*mCodonFreq[21]-aOmega*mCodonFreq[29]-aOmega*mCodonFreq[33]-aK*mCodonFreq[38]-aOmega*mCodonFreq[39]-aOmega*mCodonFreq[40]-aK*aOmega*mCodonFreq[41]-aK*aOmega*mCodonFreq[53];
	mQ[2295] = aK*mCodonFreq[38];
	mQ[2296] = aOmega*mCodonFreq[39];
	mQ[2297] = aOmega*mCodonFreq[40];
	mQ[2298] = aK*aOmega*mCodonFreq[41];
	mQ[2310] = aK*aOmega*mCodonFreq[53];
	mQ[2327] = aOmega*mCodonFreq[9];
	mQ[2340] = aOmega*mCodonFreq[22];
	mQ[2348] = aOmega*mCodonFreq[30];
	mQ[2352] = aOmega*mCodonFreq[34];
	mQ[2355] = aK*mCodonFreq[37];
	mQ[2356] = -aOmega*mCodonFreq[9]-aOmega*mCodonFreq[22]-aOmega*mCodonFreq[30]-aOmega*mCodonFreq[34]-aK*mCodonFreq[37]-aOmega*mCodonFreq[39]-aOmega*mCodonFreq[40]-aK*aOmega*mCodonFreq[42]-aK*aOmega*mCodonFreq[54];
	mQ[2357] = aOmega*mCodonFreq[39];
	mQ[2358] = aOmega*mCodonFreq[40];
	mQ[2360] = aK*aOmega*mCodonFreq[42];
	mQ[2372] = aK*aOmega*mCodonFreq[54];
	mQ[2402] = aOmega*mCodonFreq[23];
	mQ[2410] = aOmega*mCodonFreq[31];
	mQ[2414] = aOmega*mCodonFreq[35];
	mQ[2416] = aOmega*mCodonFreq[37];
	mQ[2417] = aOmega*mCodonFreq[38];
	mQ[2418] = -aOmega*mCodonFreq[23]-aOmega*mCodonFreq[31]-aOmega*mCodonFreq[35]-aOmega*mCodonFreq[37]-aOmega*mCodonFreq[38]-aK*mCodonFreq[40]-aK*aOmega*mCodonFreq[43]-aK*aOmega*mCodonFreq[55];
	mQ[2419] = aK*mCodonFreq[40];
	mQ[2422] = aK*aOmega*mCodonFreq[43];
	mQ[2434] = aK*aOmega*mCodonFreq[55];
	mQ[2464] = aOmega*mCodonFreq[24];
	mQ[2472] = aOmega*mCodonFreq[32];
	mQ[2476] = aOmega*mCodonFreq[36];
	mQ[2477] = aOmega*mCodonFreq[37];
	mQ[2478] = aOmega*mCodonFreq[38];
	mQ[2479] = aK*mCodonFreq[39];
	mQ[2480] = -aOmega*mCodonFreq[24]-aOmega*mCodonFreq[32]-aOmega*mCodonFreq[36]-aOmega*mCodonFreq[37]-aOmega*mCodonFreq[38]-aK*mCodonFreq[39]-aK*aOmega*mCodonFreq[44]-aK*aOmega*mCodonFreq[56];
	mQ[2484] = aK*aOmega*mCodonFreq[44];
	mQ[2496] = aK*aOmega*mCodonFreq[56];
	mQ[2511] = aOmega*mCodonFreq[10];
	mQ[2526] = aOmega*mCodonFreq[25];
	mQ[2530] = aOmega*mCodonFreq[29];
	mQ[2534] = aOmega*mCodonFreq[33];
	mQ[2538] = aK*aOmega*mCodonFreq[37];
	mQ[2542] = -aOmega*mCodonFreq[10]-aOmega*mCodonFreq[25]-aOmega*mCodonFreq[29]-aOmega*mCodonFreq[33]-aK*aOmega*mCodonFreq[37]-aK*mCodonFreq[42]-aOmega*mCodonFreq[43]-aOmega*mCodonFreq[44]-aK*aOmega*mCodonFreq[57];
	mQ[2543] = aK*mCodonFreq[42];
	mQ[2544] = aOmega*mCodonFreq[43];
	mQ[2545] = aOmega*mCodonFreq[44];
	mQ[2558] = aK*aOmega*mCodonFreq[57];
	mQ[2573] = aOmega*mCodonFreq[11];
	mQ[2588] = aOmega*mCodonFreq[26];
	mQ[2592] = aOmega*mCodonFreq[30];
	mQ[2596] = aOmega*mCodonFreq[34];
	mQ[2600] = aK*aOmega*mCodonFreq[38];
	mQ[2603] = aK*mCodonFreq[41];
	mQ[2604] = -aOmega*mCodonFreq[11]-aOmega*mCodonFreq[26]-aOmega*mCodonFreq[30]-aOmega*mCodonFreq[34]-aK*aOmega*mCodonFreq[38]-aK*mCodonFreq[41]-aOmega*mCodonFreq[43]-aOmega*mCodonFreq[44]-aK*aOmega*mCodonFreq[58];
	mQ[2605] = aOmega*mCodonFreq[43];
	mQ[2606] = aOmega*mCodonFreq[44];
	mQ[2620] = aK*aOmega*mCodonFreq[58];
	mQ[2650] = mCodonFreq[27];
	mQ[2654] = aOmega*mCodonFreq[31];
	mQ[2658] = aOmega*mCodonFreq[35];
	mQ[2662] = aK*aOmega*mCodonFreq[39];
	mQ[2664] = aOmega*mCodonFreq[41];
	mQ[2665] = aOmega*mCodonFreq[42];
	mQ[2666] = -mCodonFreq[27]-aOmega*mCodonFreq[31]-aOmega*mCodonFreq[35]-aK*aOmega*mCodonFreq[39]-aOmega*mCodonFreq[41]-aOmega*mCodonFreq[42]-aK*mCodonFreq[44]-aK*aOmega*mCodonFreq[59];
	mQ[2667] = aK*mCodonFreq[44];
	mQ[2682] = aK*aOmega*mCodonFreq[59];
	mQ[2696] = aOmega*mCodonFreq[12];
	mQ[2712] = mCodonFreq[28];
	mQ[2716] = aOmega*mCodonFreq[32];
	mQ[2720] = aOmega*mCodonFreq[36];
	mQ[2724] = aK*aOmega*mCodonFreq[40];
	mQ[2725] = aOmega*mCodonFreq[41];
	mQ[2726] = aOmega*mCodonFreq[42];
	mQ[2727] = aK*mCodonFreq[43];
	mQ[2728] = -aOmega*mCodonFreq[12]-mCodonFreq[28]-aOmega*mCodonFreq[32]-aOmega*mCodonFreq[36]-aK*aOmega*mCodonFreq[40]-aOmega*mCodonFreq[41]-aOmega*mCodonFreq[42]-aK*mCodonFreq[43]-aK*aOmega*mCodonFreq[60];
	mQ[2744] = aK*aOmega*mCodonFreq[60];
	mQ[2745] = aOmega*mCodonFreq[0];
	mQ[2758] = aOmega*mCodonFreq[13];
	mQ[2774] = aK*aOmega*mCodonFreq[29];
	mQ[2790] = -aOmega*mCodonFreq[0]-aOmega*mCodonFreq[13]-aK*aOmega*mCodonFreq[29]-aK*mCodonFreq[46]-mCodonFreq[47]-mCodonFreq[48]-aK*aOmega*mCodonFreq[49]-aOmega*mCodonFreq[53]-aOmega*mCodonFreq[57];
	mQ[2791] = aK*mCodonFreq[46];
	mQ[2792] = mCodonFreq[47];
	mQ[2793] = mCodonFreq[48];
	mQ[2794] = aK*aOmega*mCodonFreq[49];
	mQ[2798] = aOmega*mCodonFreq[53];
	mQ[2802] = aOmega*mCodonFreq[57];
	mQ[2807] = aOmega*mCodonFreq[1];
	mQ[2820] = aOmega*mCodonFreq[14];
	mQ[2836] = aK*aOmega*mCodonFreq[30];
	mQ[2851] = aK*mCodonFreq[45];
	mQ[2852] = -aOmega*mCodonFreq[1]-aOmega*mCodonFreq[14]-aK*aOmega*mCodonFreq[30]-aK*mCodonFreq[45]-mCodonFreq[47]-mCodonFreq[48]-aK*aOmega*mCodonFreq[50]-aOmega*mCodonFreq[54]-aOmega*mCodonFreq[58];
	mQ[2853] = mCodonFreq[47];
	mQ[2854] = mCodonFreq[48];
	mQ[2856] = aK*aOmega*mCodonFreq[50];
	mQ[2860] = aOmega*mCodonFreq[54];
	mQ[2864] = aOmega*mCodonFreq[58];
	mQ[2869] = aOmega*mCodonFreq[2];
	mQ[2882] = aOmega*mCodonFreq[15];
	mQ[2898] = aK*aOmega*mCodonFreq[31];
	mQ[2912] = mCodonFreq[45];
	mQ[2913] = mCodonFreq[46];
	mQ[2914] = -aOmega*mCodonFreq[2]-aOmega*mCodonFreq[15]-aK*aOmega*mCodonFreq[31]-mCodonFreq[45]-mCodonFreq[46]-aK*mCodonFreq[48]-aK*aOmega*mCodonFreq[51]-aOmega*mCodonFreq[55]-aOmega*mCodonFreq[59];
	mQ[2915] = aK*mCodonFreq[48];
	mQ[2918] = aK*aOmega*mCodonFreq[51];
	mQ[2922] = aOmega*mCodonFreq[55];
	mQ[2926] = aOmega*mCodonFreq[59];
	mQ[2931] = aOmega*mCodonFreq[3];
	mQ[2944] = aOmega*mCodonFreq[16];
	mQ[2960] = aK*aOmega*mCodonFreq[32];
	mQ[2973] = mCodonFreq[45];
	mQ[2974] = mCodonFreq[46];
	mQ[2975] = aK*mCodonFreq[47];
	mQ[2976] = -aOmega*mCodonFreq[3]-aOmega*mCodonFreq[16]-aK*aOmega*mCodonFreq[32]-mCodonFreq[45]-mCodonFreq[46]-aK*mCodonFreq[47]-aK*aOmega*mCodonFreq[52]-aOmega*mCodonFreq[56]-aOmega*mCodonFreq[60];
	mQ[2980] = aK*aOmega*mCodonFreq[52];
	mQ[2984] = aOmega*mCodonFreq[56];
	mQ[2988] = aOmega*mCodonFreq[60];
	mQ[2993] = aOmega*mCodonFreq[4];
	mQ[3006] = aOmega*mCodonFreq[17];
	mQ[3022] = aK*aOmega*mCodonFreq[33];
	mQ[3034] = aK*aOmega*mCodonFreq[45];
	mQ[3038] = -aOmega*mCodonFreq[4]-aOmega*mCodonFreq[17]-aK*aOmega*mCodonFreq[33]-aK*aOmega*mCodonFreq[45]-aK*mCodonFreq[50]-mCodonFreq[51]-mCodonFreq[52]-aOmega*mCodonFreq[53]-aOmega*mCodonFreq[57];
	mQ[3039] = aK*mCodonFreq[50];
	mQ[3040] = mCodonFreq[51];
	mQ[3041] = mCodonFreq[52];
	mQ[3042] = aOmega*mCodonFreq[53];
	mQ[3046] = aOmega*mCodonFreq[57];
	mQ[3055] = aOmega*mCodonFreq[5];
	mQ[3068] = aOmega*mCodonFreq[18];
	mQ[3084] = aK*aOmega*mCodonFreq[34];
	mQ[3096] = aK*aOmega*mCodonFreq[46];
	mQ[3099] = aK*mCodonFreq[49];
	mQ[3100] = -aOmega*mCodonFreq[5]-aOmega*mCodonFreq[18]-aK*aOmega*mCodonFreq[34]-aK*aOmega*mCodonFreq[46]-aK*mCodonFreq[49]-mCodonFreq[51]-mCodonFreq[52]-aOmega*mCodonFreq[54]-aOmega*mCodonFreq[58];
	mQ[3101] = mCodonFreq[51];
	mQ[3102] = mCodonFreq[52];
	mQ[3104] = aOmega*mCodonFreq[54];
	mQ[3108] = aOmega*mCodonFreq[58];
	mQ[3117] = aOmega*mCodonFreq[6];
	mQ[3130] = aOmega*mCodonFreq[19];
	mQ[3146] = aK*aOmega*mCodonFreq[35];
	mQ[3158] = aK*aOmega*mCodonFreq[47];
	mQ[3160] = mCodonFreq[49];
	mQ[3161] = mCodonFreq[50];
	mQ[3162] = -aOmega*mCodonFreq[6]-aOmega*mCodonFreq[19]-aK*aOmega*mCodonFreq[35]-aK*aOmega*mCodonFreq[47]-mCodonFreq[49]-mCodonFreq[50]-aK*mCodonFreq[52]-aOmega*mCodonFreq[55]-aOmega*mCodonFreq[59];
	mQ[3163] = aK*mCodonFreq[52];
	mQ[3166] = aOmega*mCodonFreq[55];
	mQ[3170] = aOmega*mCodonFreq[59];
	mQ[3179] = aOmega*mCodonFreq[7];
	mQ[3192] = aOmega*mCodonFreq[20];
	mQ[3208] = aK*aOmega*mCodonFreq[36];
	mQ[3220] = aK*aOmega*mCodonFreq[48];
	mQ[3221] = mCodonFreq[49];
	mQ[3222] = mCodonFreq[50];
	mQ[3223] = aK*mCodonFreq[51];
	mQ[3224] = -aOmega*mCodonFreq[7]-aOmega*mCodonFreq[20]-aK*aOmega*mCodonFreq[36]-aK*aOmega*mCodonFreq[48]-mCodonFreq[49]-mCodonFreq[50]-aK*mCodonFreq[51]-aOmega*mCodonFreq[56]-aOmega*mCodonFreq[60];
	mQ[3228] = aOmega*mCodonFreq[56];
	mQ[3232] = aOmega*mCodonFreq[60];
	mQ[3241] = aOmega*mCodonFreq[8];
	mQ[3254] = aOmega*mCodonFreq[21];
	mQ[3270] = aK*aOmega*mCodonFreq[37];
	mQ[3278] = aOmega*mCodonFreq[45];
	mQ[3282] = aOmega*mCodonFreq[49];
	mQ[3286] = -aOmega*mCodonFreq[8]-aOmega*mCodonFreq[21]-aK*aOmega*mCodonFreq[37]-aOmega*mCodonFreq[45]-aOmega*mCodonFreq[49]-aK*mCodonFreq[54]-aOmega*mCodonFreq[55]-aOmega*mCodonFreq[56]-aK*aOmega*mCodonFreq[57];
	mQ[3287] = aK*mCodonFreq[54];
	mQ[3288] = aOmega*mCodonFreq[55];
	mQ[3289] = aOmega*mCodonFreq[56];
	mQ[3290] = aK*aOmega*mCodonFreq[57];
	mQ[3303] = aOmega*mCodonFreq[9];
	mQ[3316] = aOmega*mCodonFreq[22];
	mQ[3332] = aK*aOmega*mCodonFreq[38];
	mQ[3340] = aOmega*mCodonFreq[46];
	mQ[3344] = aOmega*mCodonFreq[50];
	mQ[3347] = aK*mCodonFreq[53];
	mQ[3348] = -aOmega*mCodonFreq[9]-aOmega*mCodonFreq[22]-aK*aOmega*mCodonFreq[38]-aOmega*mCodonFreq[46]-aOmega*mCodonFreq[50]-aK*mCodonFreq[53]-aOmega*mCodonFreq[55]-aOmega*mCodonFreq[56]-aK*aOmega*mCodonFreq[58];
	mQ[3349] = aOmega*mCodonFreq[55];
	mQ[3350] = aOmega*mCodonFreq[56];
	mQ[3352] = aK*aOmega*mCodonFreq[58];
	mQ[3378] = aOmega*mCodonFreq[23];
	mQ[3394] = aK*aOmega*mCodonFreq[39];
	mQ[3402] = aOmega*mCodonFreq[47];
	mQ[3406] = aOmega*mCodonFreq[51];
	mQ[3408] = aOmega*mCodonFreq[53];
	mQ[3409] = aOmega*mCodonFreq[54];
	mQ[3410] = -aOmega*mCodonFreq[23]-aK*aOmega*mCodonFreq[39]-aOmega*mCodonFreq[47]-aOmega*mCodonFreq[51]-aOmega*mCodonFreq[53]-aOmega*mCodonFreq[54]-aK*mCodonFreq[56]-aK*aOmega*mCodonFreq[59];
	mQ[3411] = aK*mCodonFreq[56];
	mQ[3414] = aK*aOmega*mCodonFreq[59];
	mQ[3440] = aOmega*mCodonFreq[24];
	mQ[3456] = aK*aOmega*mCodonFreq[40];
	mQ[3464] = aOmega*mCodonFreq[48];
	mQ[3468] = aOmega*mCodonFreq[52];
	mQ[3469] = aOmega*mCodonFreq[53];
	mQ[3470] = aOmega*mCodonFreq[54];
	mQ[3471] = aK*mCodonFreq[55];
	mQ[3472] = -aOmega*mCodonFreq[24]-aK*aOmega*mCodonFreq[40]-aOmega*mCodonFreq[48]-aOmega*mCodonFreq[52]-aOmega*mCodonFreq[53]-aOmega*mCodonFreq[54]-aK*mCodonFreq[55]-aK*aOmega*mCodonFreq[60];
	mQ[3476] = aK*aOmega*mCodonFreq[60];
	mQ[3487] = aOmega*mCodonFreq[10];
	mQ[3502] = aOmega*mCodonFreq[25];
	mQ[3518] = aK*aOmega*mCodonFreq[41];
	mQ[3522] = aOmega*mCodonFreq[45];
	mQ[3526] = aOmega*mCodonFreq[49];
	mQ[3530] = aK*aOmega*mCodonFreq[53];
	mQ[3534] = -aOmega*mCodonFreq[10]-aOmega*mCodonFreq[25]-aK*aOmega*mCodonFreq[41]-aOmega*mCodonFreq[45]-aOmega*mCodonFreq[49]-aK*aOmega*mCodonFreq[53]-aK*mCodonFreq[58]-mCodonFreq[59]-mCodonFreq[60];
	mQ[3535] = aK*mCodonFreq[58];
	mQ[3536] = mCodonFreq[59];
	mQ[3537] = mCodonFreq[60];
	mQ[3549] = aOmega*mCodonFreq[11];
	mQ[3564] = aOmega*mCodonFreq[26];
	mQ[3580] = aK*aOmega*mCodonFreq[42];
	mQ[3584] = aOmega*mCodonFreq[46];
	mQ[3588] = aOmega*mCodonFreq[50];
	mQ[3592] = aK*aOmega*mCodonFreq[54];
	mQ[3595] = aK*mCodonFreq[57];
	mQ[3596] = -aOmega*mCodonFreq[11]-aOmega*mCodonFreq[26]-aK*aOmega*mCodonFreq[42]-aOmega*mCodonFreq[46]-aOmega*mCodonFreq[50]-aK*aOmega*mCodonFreq[54]-aK*mCodonFreq[57]-mCodonFreq[59]-mCodonFreq[60];
	mQ[3597] = mCodonFreq[59];
	mQ[3598] = mCodonFreq[60];
	mQ[3626] = aOmega*mCodonFreq[27];
	mQ[3642] = aK*aOmega*mCodonFreq[43];
	mQ[3646] = aOmega*mCodonFreq[47];
	mQ[3650] = aOmega*mCodonFreq[51];
	mQ[3654] = aK*aOmega*mCodonFreq[55];
	mQ[3656] = mCodonFreq[57];
	mQ[3657] = mCodonFreq[58];
	mQ[3658] = -aOmega*mCodonFreq[27]-aK*aOmega*mCodonFreq[43]-aOmega*mCodonFreq[47]-aOmega*mCodonFreq[51]-aK*aOmega*mCodonFreq[55]-mCodonFreq[57]-mCodonFreq[58]-aK*mCodonFreq[60];
	mQ[3659] = aK*mCodonFreq[60];
	mQ[3672] = aOmega*mCodonFreq[12];
	mQ[3688] = aOmega*mCodonFreq[28];
	mQ[3704] = aK*aOmega*mCodonFreq[44];
	mQ[3708] = aOmega*mCodonFreq[48];
	mQ[3712] = aOmega*mCodonFreq[52];
	mQ[3716] = aK*aOmega*mCodonFreq[56];
	mQ[3717] = mCodonFreq[57];
	mQ[3718] = mCodonFreq[58];
	mQ[3719] = aK*mCodonFreq[59];
	mQ[3720] = -aOmega*mCodonFreq[12]-aOmega*mCodonFreq[28]-aK*aOmega*mCodonFreq[44]-aOmega*mCodonFreq[48]-aOmega*mCodonFreq[52]-aK*aOmega*mCodonFreq[56]-mCodonFreq[57]-mCodonFreq[58]-aK*mCodonFreq[59];

	// Compute the scale factor
	double scale_q =     aK*mCodonFreq[0]*mCodonFreq[1]
	             + aOmega*mCodonFreq[0]*mCodonFreq[2]
	             + aOmega*mCodonFreq[0]*mCodonFreq[3]
	             + aK*aOmega*mCodonFreq[0]*mCodonFreq[4]
	             + aOmega*mCodonFreq[0]*mCodonFreq[8]
	             + aOmega*mCodonFreq[0]*mCodonFreq[10]
	             + aK*aOmega*mCodonFreq[0]*mCodonFreq[13]
	             + aOmega*mCodonFreq[0]*mCodonFreq[29]
	             + aOmega*mCodonFreq[0]*mCodonFreq[45]
	             + aOmega*mCodonFreq[1]*mCodonFreq[2]
	             + aOmega*mCodonFreq[1]*mCodonFreq[3]
	             + aK*aOmega*mCodonFreq[1]*mCodonFreq[5]
	             + aOmega*mCodonFreq[1]*mCodonFreq[9]
	             + aOmega*mCodonFreq[1]*mCodonFreq[11]
	             + aK*aOmega*mCodonFreq[1]*mCodonFreq[14]
	             + aOmega*mCodonFreq[1]*mCodonFreq[30]
	             + aOmega*mCodonFreq[1]*mCodonFreq[46]
	             + aK*mCodonFreq[2]*mCodonFreq[3]
	             + aK*aOmega*mCodonFreq[2]*mCodonFreq[6]
	             + aK*mCodonFreq[2]*mCodonFreq[15]
	             + aOmega*mCodonFreq[2]*mCodonFreq[31]
	             + aOmega*mCodonFreq[2]*mCodonFreq[47]
	             + aK*aOmega*mCodonFreq[3]*mCodonFreq[7]
	             + aOmega*mCodonFreq[3]*mCodonFreq[12]
	             + aK*mCodonFreq[3]*mCodonFreq[16]
	             + aOmega*mCodonFreq[3]*mCodonFreq[32]
	             + aOmega*mCodonFreq[3]*mCodonFreq[48]
	             + aK*mCodonFreq[4]*mCodonFreq[5]
	             +   mCodonFreq[4]*mCodonFreq[6]
	             +   mCodonFreq[4]*mCodonFreq[7]
	             + aOmega*mCodonFreq[4]*mCodonFreq[8]
	             + aOmega*mCodonFreq[4]*mCodonFreq[10]
	             + aK*aOmega*mCodonFreq[4]*mCodonFreq[17]
	             + aOmega*mCodonFreq[4]*mCodonFreq[33]
	             + aOmega*mCodonFreq[4]*mCodonFreq[49]
	             +   mCodonFreq[5]*mCodonFreq[6]
	             +   mCodonFreq[5]*mCodonFreq[7]
	             + aOmega*mCodonFreq[5]*mCodonFreq[9]
	             + aOmega*mCodonFreq[5]*mCodonFreq[11]
	             + aK*aOmega*mCodonFreq[5]*mCodonFreq[18]
	             + aOmega*mCodonFreq[5]*mCodonFreq[34]
	             + aOmega*mCodonFreq[5]*mCodonFreq[50]
	             + aK*mCodonFreq[6]*mCodonFreq[7]
	             + aK*aOmega*mCodonFreq[6]*mCodonFreq[19]
	             + aOmega*mCodonFreq[6]*mCodonFreq[35]
	             + aOmega*mCodonFreq[6]*mCodonFreq[51]
	             + aOmega*mCodonFreq[7]*mCodonFreq[12]
	             + aK*aOmega*mCodonFreq[7]*mCodonFreq[20]
	             + aOmega*mCodonFreq[7]*mCodonFreq[36]
	             + aOmega*mCodonFreq[7]*mCodonFreq[52]
	             + aK*mCodonFreq[8]*mCodonFreq[9]
	             + aK*aOmega*mCodonFreq[8]*mCodonFreq[10]
	             + aK*aOmega*mCodonFreq[8]*mCodonFreq[21]
	             + aOmega*mCodonFreq[8]*mCodonFreq[37]
	             + aOmega*mCodonFreq[8]*mCodonFreq[53]
	             + aK*aOmega*mCodonFreq[9]*mCodonFreq[11]
	             + aK*aOmega*mCodonFreq[9]*mCodonFreq[22]
	             + aOmega*mCodonFreq[9]*mCodonFreq[38]
	             + aOmega*mCodonFreq[9]*mCodonFreq[54]
	             + aK*mCodonFreq[10]*mCodonFreq[11]
	             + aOmega*mCodonFreq[10]*mCodonFreq[12]
	             + aK*aOmega*mCodonFreq[10]*mCodonFreq[25]
	             + aOmega*mCodonFreq[10]*mCodonFreq[41]
	             + aOmega*mCodonFreq[10]*mCodonFreq[57]
	             + aOmega*mCodonFreq[11]*mCodonFreq[12]
	             + aK*aOmega*mCodonFreq[11]*mCodonFreq[26]
	             + aOmega*mCodonFreq[11]*mCodonFreq[42]
	             + aOmega*mCodonFreq[11]*mCodonFreq[58]
	             + aK*aOmega*mCodonFreq[12]*mCodonFreq[28]
	             + aOmega*mCodonFreq[12]*mCodonFreq[44]
	             + aOmega*mCodonFreq[12]*mCodonFreq[60]
	             + aK*mCodonFreq[13]*mCodonFreq[14]
	             +   mCodonFreq[13]*mCodonFreq[15]
	             +   mCodonFreq[13]*mCodonFreq[16]
	             + aK*aOmega*mCodonFreq[13]*mCodonFreq[17]
	             + aOmega*mCodonFreq[13]*mCodonFreq[21]
	             + aOmega*mCodonFreq[13]*mCodonFreq[25]
	             + aOmega*mCodonFreq[13]*mCodonFreq[29]
	             + aOmega*mCodonFreq[13]*mCodonFreq[45]
	             +   mCodonFreq[14]*mCodonFreq[15]
	             +   mCodonFreq[14]*mCodonFreq[16]
	             + aK*aOmega*mCodonFreq[14]*mCodonFreq[18]
	             + aOmega*mCodonFreq[14]*mCodonFreq[22]
	             + aOmega*mCodonFreq[14]*mCodonFreq[26]
	             + aOmega*mCodonFreq[14]*mCodonFreq[30]
	             + aOmega*mCodonFreq[14]*mCodonFreq[46]
	             + aK*mCodonFreq[15]*mCodonFreq[16]
	             + aK*aOmega*mCodonFreq[15]*mCodonFreq[19]
	             + aOmega*mCodonFreq[15]*mCodonFreq[23]
	             + aOmega*mCodonFreq[15]*mCodonFreq[27]
	             + aOmega*mCodonFreq[15]*mCodonFreq[31]
	             + aOmega*mCodonFreq[15]*mCodonFreq[47]
	             + aK*aOmega*mCodonFreq[16]*mCodonFreq[20]
	             + aOmega*mCodonFreq[16]*mCodonFreq[24]
	             + aOmega*mCodonFreq[16]*mCodonFreq[28]
	             + aOmega*mCodonFreq[16]*mCodonFreq[32]
	             + aOmega*mCodonFreq[16]*mCodonFreq[48]
	             + aK*mCodonFreq[17]*mCodonFreq[18]
	             +   mCodonFreq[17]*mCodonFreq[19]
	             +   mCodonFreq[17]*mCodonFreq[20]
	             + aOmega*mCodonFreq[17]*mCodonFreq[21]
	             + aOmega*mCodonFreq[17]*mCodonFreq[25]
	             + aOmega*mCodonFreq[17]*mCodonFreq[33]
	             + aOmega*mCodonFreq[17]*mCodonFreq[49]
	             +   mCodonFreq[18]*mCodonFreq[19]
	             +   mCodonFreq[18]*mCodonFreq[20]
	             + aOmega*mCodonFreq[18]*mCodonFreq[22]
	             + aOmega*mCodonFreq[18]*mCodonFreq[26]
	             + aOmega*mCodonFreq[18]*mCodonFreq[34]
	             + aOmega*mCodonFreq[18]*mCodonFreq[50]
	             + aK*mCodonFreq[19]*mCodonFreq[20]
	             + aOmega*mCodonFreq[19]*mCodonFreq[23]
	             + aOmega*mCodonFreq[19]*mCodonFreq[27]
	             + aOmega*mCodonFreq[19]*mCodonFreq[35]
	             + aOmega*mCodonFreq[19]*mCodonFreq[51]
	             + aOmega*mCodonFreq[20]*mCodonFreq[24]
	             + aOmega*mCodonFreq[20]*mCodonFreq[28]
	             + aOmega*mCodonFreq[20]*mCodonFreq[36]
	             + aOmega*mCodonFreq[20]*mCodonFreq[52]
	             + aK*mCodonFreq[21]*mCodonFreq[22]
	             + aOmega*mCodonFreq[21]*mCodonFreq[23]
	             + aOmega*mCodonFreq[21]*mCodonFreq[24]
	             + aK*aOmega*mCodonFreq[21]*mCodonFreq[25]
	             + aOmega*mCodonFreq[21]*mCodonFreq[37]
	             + aOmega*mCodonFreq[21]*mCodonFreq[53]
	             + aOmega*mCodonFreq[22]*mCodonFreq[23]
	             + aOmega*mCodonFreq[22]*mCodonFreq[24]
	             + aK*aOmega*mCodonFreq[22]*mCodonFreq[26]
	             + aOmega*mCodonFreq[22]*mCodonFreq[38]
	             + aOmega*mCodonFreq[22]*mCodonFreq[54]
	             + aK*mCodonFreq[23]*mCodonFreq[24]
	             + aK*aOmega*mCodonFreq[23]*mCodonFreq[27]
	             + aOmega*mCodonFreq[23]*mCodonFreq[39]
	             + aOmega*mCodonFreq[23]*mCodonFreq[55]
	             + aK*aOmega*mCodonFreq[24]*mCodonFreq[28]
	             + aOmega*mCodonFreq[24]*mCodonFreq[40]
	             + aOmega*mCodonFreq[24]*mCodonFreq[56]
	             + aK*mCodonFreq[25]*mCodonFreq[26]
	             +   mCodonFreq[25]*mCodonFreq[27]
	             +   mCodonFreq[25]*mCodonFreq[28]
	             + aOmega*mCodonFreq[25]*mCodonFreq[41]
	             + aOmega*mCodonFreq[25]*mCodonFreq[57]
	             +   mCodonFreq[26]*mCodonFreq[27]
	             +   mCodonFreq[26]*mCodonFreq[28]
	             + aOmega*mCodonFreq[26]*mCodonFreq[42]
	             + aOmega*mCodonFreq[26]*mCodonFreq[58]
	             + aK*mCodonFreq[27]*mCodonFreq[28]
	             +   mCodonFreq[27]*mCodonFreq[43]
	             + aOmega*mCodonFreq[27]*mCodonFreq[59]
	             +   mCodonFreq[28]*mCodonFreq[44]
	             + aOmega*mCodonFreq[28]*mCodonFreq[60]
	             + aK*mCodonFreq[29]*mCodonFreq[30]
	             +   mCodonFreq[29]*mCodonFreq[31]
	             + aOmega*mCodonFreq[29]*mCodonFreq[32]
	             + aK*aOmega*mCodonFreq[29]*mCodonFreq[33]
	             + aOmega*mCodonFreq[29]*mCodonFreq[37]
	             + aOmega*mCodonFreq[29]*mCodonFreq[41]
	             + aK*aOmega*mCodonFreq[29]*mCodonFreq[45]
	             +   mCodonFreq[30]*mCodonFreq[31]
	             + aOmega*mCodonFreq[30]*mCodonFreq[32]
	             + aK*aOmega*mCodonFreq[30]*mCodonFreq[34]
	             + aOmega*mCodonFreq[30]*mCodonFreq[38]
	             + aOmega*mCodonFreq[30]*mCodonFreq[42]
	             + aK*aOmega*mCodonFreq[30]*mCodonFreq[46]
	             + aK*aOmega*mCodonFreq[31]*mCodonFreq[32]
	             + aK*aOmega*mCodonFreq[31]*mCodonFreq[35]
	             + aOmega*mCodonFreq[31]*mCodonFreq[39]
	             + aOmega*mCodonFreq[31]*mCodonFreq[43]
	             + aK*aOmega*mCodonFreq[31]*mCodonFreq[47]
	             + aK*aOmega*mCodonFreq[32]*mCodonFreq[36]
	             + aOmega*mCodonFreq[32]*mCodonFreq[40]
	             + aOmega*mCodonFreq[32]*mCodonFreq[44]
	             + aK*aOmega*mCodonFreq[32]*mCodonFreq[48]
	             + aK*mCodonFreq[33]*mCodonFreq[34]
	             +   mCodonFreq[33]*mCodonFreq[35]
	             +   mCodonFreq[33]*mCodonFreq[36]
	             + aOmega*mCodonFreq[33]*mCodonFreq[37]
	             + aOmega*mCodonFreq[33]*mCodonFreq[41]
	             + aK*aOmega*mCodonFreq[33]*mCodonFreq[49]
	             +   mCodonFreq[34]*mCodonFreq[35]
	             +   mCodonFreq[34]*mCodonFreq[36]
	             + aOmega*mCodonFreq[34]*mCodonFreq[38]
	             + aOmega*mCodonFreq[34]*mCodonFreq[42]
	             + aK*aOmega*mCodonFreq[34]*mCodonFreq[50]
	             + aK*mCodonFreq[35]*mCodonFreq[36]
	             + aOmega*mCodonFreq[35]*mCodonFreq[39]
	             + aOmega*mCodonFreq[35]*mCodonFreq[43]
	             + aK*aOmega*mCodonFreq[35]*mCodonFreq[51]
	             + aOmega*mCodonFreq[36]*mCodonFreq[40]
	             + aOmega*mCodonFreq[36]*mCodonFreq[44]
	             + aK*aOmega*mCodonFreq[36]*mCodonFreq[52]
	             + aK*mCodonFreq[37]*mCodonFreq[38]
	             + aOmega*mCodonFreq[37]*mCodonFreq[39]
	             + aOmega*mCodonFreq[37]*mCodonFreq[40]
	             + aK*aOmega*mCodonFreq[37]*mCodonFreq[41]
	             + aK*aOmega*mCodonFreq[37]*mCodonFreq[53]
	             + aOmega*mCodonFreq[38]*mCodonFreq[39]
	             + aOmega*mCodonFreq[38]*mCodonFreq[40]
	             + aK*aOmega*mCodonFreq[38]*mCodonFreq[42]
	             + aK*aOmega*mCodonFreq[38]*mCodonFreq[54]
	             + aK*mCodonFreq[39]*mCodonFreq[40]
	             + aK*aOmega*mCodonFreq[39]*mCodonFreq[43]
	             + aK*aOmega*mCodonFreq[39]*mCodonFreq[55]
	             + aK*aOmega*mCodonFreq[40]*mCodonFreq[44]
	             + aK*aOmega*mCodonFreq[40]*mCodonFreq[56]
	             + aK*mCodonFreq[41]*mCodonFreq[42]
	             + aOmega*mCodonFreq[41]*mCodonFreq[43]
	             + aOmega*mCodonFreq[41]*mCodonFreq[44]
	             + aK*aOmega*mCodonFreq[41]*mCodonFreq[57]
	             + aOmega*mCodonFreq[42]*mCodonFreq[43]
	             + aOmega*mCodonFreq[42]*mCodonFreq[44]
	             + aK*aOmega*mCodonFreq[42]*mCodonFreq[58]
	             + aK*mCodonFreq[43]*mCodonFreq[44]
	             + aK*aOmega*mCodonFreq[43]*mCodonFreq[59]
	             + aK*aOmega*mCodonFreq[44]*mCodonFreq[60]
	             + aK*mCodonFreq[45]*mCodonFreq[46]
	             +   mCodonFreq[45]*mCodonFreq[47]
	             +   mCodonFreq[45]*mCodonFreq[48]
	             + aK*aOmega*mCodonFreq[45]*mCodonFreq[49]
	             + aOmega*mCodonFreq[45]*mCodonFreq[53]
	             + aOmega*mCodonFreq[45]*mCodonFreq[57]
	             +   mCodonFreq[46]*mCodonFreq[47]
	             +   mCodonFreq[46]*mCodonFreq[48]
	             + aK*aOmega*mCodonFreq[46]*mCodonFreq[50]
	             + aOmega*mCodonFreq[46]*mCodonFreq[54]
	             + aOmega*mCodonFreq[46]*mCodonFreq[58]
	             + aK*mCodonFreq[47]*mCodonFreq[48]
	             + aK*aOmega*mCodonFreq[47]*mCodonFreq[51]
	             + aOmega*mCodonFreq[47]*mCodonFreq[55]
	             + aOmega*mCodonFreq[47]*mCodonFreq[59]
	             + aK*aOmega*mCodonFreq[48]*mCodonFreq[52]
	             + aOmega*mCodonFreq[48]*mCodonFreq[56]
	             + aOmega*mCodonFreq[48]*mCodonFreq[60]
	             + aK*mCodonFreq[49]*mCodonFreq[50]
	             +   mCodonFreq[49]*mCodonFreq[51]
	             +   mCodonFreq[49]*mCodonFreq[52]
	             + aOmega*mCodonFreq[49]*mCodonFreq[53]
	             + aOmega*mCodonFreq[49]*mCodonFreq[57]
	             +   mCodonFreq[50]*mCodonFreq[51]
	             +   mCodonFreq[50]*mCodonFreq[52]
	             + aOmega*mCodonFreq[50]*mCodonFreq[54]
	             + aOmega*mCodonFreq[50]*mCodonFreq[58]
	             + aK*mCodonFreq[51]*mCodonFreq[52]
	             + aOmega*mCodonFreq[51]*mCodonFreq[55]
	             + aOmega*mCodonFreq[51]*mCodonFreq[59]
	             + aOmega*mCodonFreq[52]*mCodonFreq[56]
	             + aOmega*mCodonFreq[52]*mCodonFreq[60]
	             + aK*mCodonFreq[53]*mCodonFreq[54]
	             + aOmega*mCodonFreq[53]*mCodonFreq[55]
	             + aOmega*mCodonFreq[53]*mCodonFreq[56]
	             + aK*aOmega*mCodonFreq[53]*mCodonFreq[57]
	             + aOmega*mCodonFreq[54]*mCodonFreq[55]
	             + aOmega*mCodonFreq[54]*mCodonFreq[56]
	             + aK*aOmega*mCodonFreq[54]*mCodonFreq[58]
	             + aK*mCodonFreq[55]*mCodonFreq[56]
	             + aK*aOmega*mCodonFreq[55]*mCodonFreq[59]
	             + aK*aOmega*mCodonFreq[56]*mCodonFreq[60]
	             + aK*mCodonFreq[57]*mCodonFreq[58]
	             +   mCodonFreq[57]*mCodonFreq[59]
	             +   mCodonFreq[57]*mCodonFreq[60]
	             +   mCodonFreq[58]*mCodonFreq[59]
	             +   mCodonFreq[58]*mCodonFreq[60]
	             + aK*mCodonFreq[59]*mCodonFreq[60];

	return scale_q*2.;
}

double TransitionMatrix::fillQ(double aK)
{
	mQ[   0] = -aK*mCodonFreq[1]-mCodonFreq[2]-mCodonFreq[3]-aK*mCodonFreq[4]-mCodonFreq[8]-mCodonFreq[10]-aK*mCodonFreq[13]-mCodonFreq[29]-mCodonFreq[45];
	mQ[   1] = aK*mCodonFreq[1];
	mQ[   2] = mCodonFreq[2];
	mQ[   3] = mCodonFreq[3];
	mQ[   4] = aK*mCodonFreq[4];
	mQ[   8] = mCodonFreq[8];
	mQ[  10] = mCodonFreq[10];
	mQ[  13] = aK*mCodonFreq[13];
	mQ[  29] = mCodonFreq[29];
	mQ[  45] = mCodonFreq[45];
	mQ[  61] = aK*mCodonFreq[0];
	mQ[  62] = -aK*mCodonFreq[0]-mCodonFreq[2]-mCodonFreq[3]-aK*mCodonFreq[5]-mCodonFreq[9]-mCodonFreq[11]-aK*mCodonFreq[14]-mCodonFreq[30]-mCodonFreq[46];
	mQ[  63] = mCodonFreq[2];
	mQ[  64] = mCodonFreq[3];
	mQ[  66] = aK*mCodonFreq[5];
	mQ[  70] = mCodonFreq[9];
	mQ[  72] = mCodonFreq[11];
	mQ[  75] = aK*mCodonFreq[14];
	mQ[  91] = mCodonFreq[30];
	mQ[ 107] = mCodonFreq[46];
	mQ[ 122] = mCodonFreq[0];
	mQ[ 123] = mCodonFreq[1];
	mQ[ 124] = -mCodonFreq[0]-mCodonFreq[1]-aK*mCodonFreq[3]-aK*mCodonFreq[6]-aK*mCodonFreq[15]-mCodonFreq[31]-mCodonFreq[47];
	mQ[ 125] = aK*mCodonFreq[3];
	mQ[ 128] = aK*mCodonFreq[6];
	mQ[ 137] = aK*mCodonFreq[15];
	mQ[ 153] = mCodonFreq[31];
	mQ[ 169] = mCodonFreq[47];
	mQ[ 183] = mCodonFreq[0];
	mQ[ 184] = mCodonFreq[1];
	mQ[ 185] = aK*mCodonFreq[2];
	mQ[ 186] = -mCodonFreq[0]-mCodonFreq[1]-aK*mCodonFreq[2]-aK*mCodonFreq[7]-mCodonFreq[12]-aK*mCodonFreq[16]-mCodonFreq[32]-mCodonFreq[48];
	mQ[ 190] = aK*mCodonFreq[7];
	mQ[ 195] = mCodonFreq[12];
	mQ[ 199] = aK*mCodonFreq[16];
	mQ[ 215] = mCodonFreq[32];
	mQ[ 231] = mCodonFreq[48];
	mQ[ 244] = aK*mCodonFreq[0];
	mQ[ 248] = -aK*mCodonFreq[0]-aK*mCodonFreq[5]-mCodonFreq[6]-mCodonFreq[7]-mCodonFreq[8]-mCodonFreq[10]-aK*mCodonFreq[17]-mCodonFreq[33]-mCodonFreq[49];
	mQ[ 249] = aK*mCodonFreq[5];
	mQ[ 250] = mCodonFreq[6];
	mQ[ 251] = mCodonFreq[7];
	mQ[ 252] = mCodonFreq[8];
	mQ[ 254] = mCodonFreq[10];
	mQ[ 261] = aK*mCodonFreq[17];
	mQ[ 277] = mCodonFreq[33];
	mQ[ 293] = mCodonFreq[49];
	mQ[ 306] = aK*mCodonFreq[1];
	mQ[ 309] = aK*mCodonFreq[4];
	mQ[ 310] = -aK*mCodonFreq[1]-aK*mCodonFreq[4]-mCodonFreq[6]-mCodonFreq[7]-mCodonFreq[9]-mCodonFreq[11]-aK*mCodonFreq[18]-mCodonFreq[34]-mCodonFreq[50];
	mQ[ 311] = mCodonFreq[6];
	mQ[ 312] = mCodonFreq[7];
	mQ[ 314] = mCodonFreq[9];
	mQ[ 316] = mCodonFreq[11];
	mQ[ 323] = aK*mCodonFreq[18];
	mQ[ 339] = mCodonFreq[34];
	mQ[ 355] = mCodonFreq[50];
	mQ[ 368] = aK*mCodonFreq[2];
	mQ[ 370] = mCodonFreq[4];
	mQ[ 371] = mCodonFreq[5];
	mQ[ 372] = -aK*mCodonFreq[2]-mCodonFreq[4]-mCodonFreq[5]-aK*mCodonFreq[7]-aK*mCodonFreq[19]-mCodonFreq[35]-mCodonFreq[51];
	mQ[ 373] = aK*mCodonFreq[7];
	mQ[ 385] = aK*mCodonFreq[19];
	mQ[ 401] = mCodonFreq[35];
	mQ[ 417] = mCodonFreq[51];
	mQ[ 430] = aK*mCodonFreq[3];
	mQ[ 431] = mCodonFreq[4];
	mQ[ 432] = mCodonFreq[5];
	mQ[ 433] = aK*mCodonFreq[6];
	mQ[ 434] = -aK*mCodonFreq[3]-mCodonFreq[4]-mCodonFreq[5]-aK*mCodonFreq[6]-mCodonFreq[12]-aK*mCodonFreq[20]-mCodonFreq[36]-mCodonFreq[52];
	mQ[ 439] = mCodonFreq[12];
	mQ[ 447] = aK*mCodonFreq[20];
	mQ[ 463] = mCodonFreq[36];
	mQ[ 479] = mCodonFreq[52];
	mQ[ 488] = mCodonFreq[0];
	mQ[ 492] = mCodonFreq[4];
	mQ[ 496] = -mCodonFreq[0]-mCodonFreq[4]-aK*mCodonFreq[9]-aK*mCodonFreq[10]-aK*mCodonFreq[21]-mCodonFreq[37]-mCodonFreq[53];
	mQ[ 497] = aK*mCodonFreq[9];
	mQ[ 498] = aK*mCodonFreq[10];
	mQ[ 509] = aK*mCodonFreq[21];
	mQ[ 525] = mCodonFreq[37];
	mQ[ 541] = mCodonFreq[53];
	mQ[ 550] = mCodonFreq[1];
	mQ[ 554] = mCodonFreq[5];
	mQ[ 557] = aK*mCodonFreq[8];
	mQ[ 558] = -mCodonFreq[1]-mCodonFreq[5]-aK*mCodonFreq[8]-aK*mCodonFreq[11]-aK*mCodonFreq[22]-mCodonFreq[38]-mCodonFreq[54];
	mQ[ 560] = aK*mCodonFreq[11];
	mQ[ 571] = aK*mCodonFreq[22];
	mQ[ 587] = mCodonFreq[38];
	mQ[ 603] = mCodonFreq[54];
	mQ[ 610] = mCodonFreq[0];
	mQ[ 614] = mCodonFreq[4];
	mQ[ 618] = aK*mCodonFreq[8];
	mQ[ 620] = -mCodonFreq[0]-mCodonFreq[4]-aK*mCodonFreq[8]-aK*mCodonFreq[11]-mCodonFreq[12]-aK*mCodonFreq[25]-mCodonFreq[41]-mCodonFreq[57];
	mQ[ 621] = aK*mCodonFreq[11];
	mQ[ 622] = mCodonFreq[12];
	mQ[ 635] = aK*mCodonFreq[25];
	mQ[ 651] = mCodonFreq[41];
	mQ[ 667] = mCodonFreq[57];
	mQ[ 672] = mCodonFreq[1];
	mQ[ 676] = mCodonFreq[5];
	mQ[ 680] = aK*mCodonFreq[9];
	mQ[ 681] = aK*mCodonFreq[10];
	mQ[ 682] = -mCodonFreq[1]-mCodonFreq[5]-aK*mCodonFreq[9]-aK*mCodonFreq[10]-mCodonFreq[12]-aK*mCodonFreq[26]-mCodonFreq[42]-mCodonFreq[58];
	mQ[ 683] = mCodonFreq[12];
	mQ[ 697] = aK*mCodonFreq[26];
	mQ[ 713] = mCodonFreq[42];
	mQ[ 729] = mCodonFreq[58];
	mQ[ 735] = mCodonFreq[3];
	mQ[ 739] = mCodonFreq[7];
	mQ[ 742] = mCodonFreq[10];
	mQ[ 743] = mCodonFreq[11];
	mQ[ 744] = -mCodonFreq[3]-mCodonFreq[7]-mCodonFreq[10]-mCodonFreq[11]-aK*mCodonFreq[28]-mCodonFreq[44]-mCodonFreq[60];
	mQ[ 760] = aK*mCodonFreq[28];
	mQ[ 776] = mCodonFreq[44];
	mQ[ 792] = mCodonFreq[60];
	mQ[ 793] = aK*mCodonFreq[0];
	mQ[ 806] = -aK*mCodonFreq[0]-aK*mCodonFreq[14]-mCodonFreq[15]-mCodonFreq[16]-aK*mCodonFreq[17]-mCodonFreq[21]-mCodonFreq[25]-mCodonFreq[29]-mCodonFreq[45];
	mQ[ 807] = aK*mCodonFreq[14];
	mQ[ 808] = mCodonFreq[15];
	mQ[ 809] = mCodonFreq[16];
	mQ[ 810] = aK*mCodonFreq[17];
	mQ[ 814] = mCodonFreq[21];
	mQ[ 818] = mCodonFreq[25];
	mQ[ 822] = mCodonFreq[29];
	mQ[ 838] = mCodonFreq[45];
	mQ[ 855] = aK*mCodonFreq[1];
	mQ[ 867] = aK*mCodonFreq[13];
	mQ[ 868] = -aK*mCodonFreq[1]-aK*mCodonFreq[13]-mCodonFreq[15]-mCodonFreq[16]-aK*mCodonFreq[18]-mCodonFreq[22]-mCodonFreq[26]-mCodonFreq[30]-mCodonFreq[46];
	mQ[ 869] = mCodonFreq[15];
	mQ[ 870] = mCodonFreq[16];
	mQ[ 872] = aK*mCodonFreq[18];
	mQ[ 876] = mCodonFreq[22];
	mQ[ 880] = mCodonFreq[26];
	mQ[ 884] = mCodonFreq[30];
	mQ[ 900] = mCodonFreq[46];
	mQ[ 917] = aK*mCodonFreq[2];
	mQ[ 928] = mCodonFreq[13];
	mQ[ 929] = mCodonFreq[14];
	mQ[ 930] = -aK*mCodonFreq[2]-mCodonFreq[13]-mCodonFreq[14]-aK*mCodonFreq[16]-aK*mCodonFreq[19]-mCodonFreq[23]-mCodonFreq[27]-mCodonFreq[31]-mCodonFreq[47];
	mQ[ 931] = aK*mCodonFreq[16];
	mQ[ 934] = aK*mCodonFreq[19];
	mQ[ 938] = mCodonFreq[23];
	mQ[ 942] = mCodonFreq[27];
	mQ[ 946] = mCodonFreq[31];
	mQ[ 962] = mCodonFreq[47];
	mQ[ 979] = aK*mCodonFreq[3];
	mQ[ 989] = mCodonFreq[13];
	mQ[ 990] = mCodonFreq[14];
	mQ[ 991] = aK*mCodonFreq[15];
	mQ[ 992] = -aK*mCodonFreq[3]-mCodonFreq[13]-mCodonFreq[14]-aK*mCodonFreq[15]-aK*mCodonFreq[20]-mCodonFreq[24]-mCodonFreq[28]-mCodonFreq[32]-mCodonFreq[48];
	mQ[ 996] = aK*mCodonFreq[20];
	mQ[1000] = mCodonFreq[24];
	mQ[1004] = mCodonFreq[28];
	mQ[1008] = mCodonFreq[32];
	mQ[1024] = mCodonFreq[48];
	mQ[1041] = aK*mCodonFreq[4];
	mQ[1050] = aK*mCodonFreq[13];
	mQ[1054] = -aK*mCodonFreq[4]-aK*mCodonFreq[13]-aK*mCodonFreq[18]-mCodonFreq[19]-mCodonFreq[20]-mCodonFreq[21]-mCodonFreq[25]-mCodonFreq[33]-mCodonFreq[49];
	mQ[1055] = aK*mCodonFreq[18];
	mQ[1056] = mCodonFreq[19];
	mQ[1057] = mCodonFreq[20];
	mQ[1058] = mCodonFreq[21];
	mQ[1062] = mCodonFreq[25];
	mQ[1070] = mCodonFreq[33];
	mQ[1086] = mCodonFreq[49];
	mQ[1103] = aK*mCodonFreq[5];
	mQ[1112] = aK*mCodonFreq[14];
	mQ[1115] = aK*mCodonFreq[17];
	mQ[1116] = -aK*mCodonFreq[5]-aK*mCodonFreq[14]-aK*mCodonFreq[17]-mCodonFreq[19]-mCodonFreq[20]-mCodonFreq[22]-mCodonFreq[26]-mCodonFreq[34]-mCodonFreq[50];
	mQ[1117] = mCodonFreq[19];
	mQ[1118] = mCodonFreq[20];
	mQ[1120] = mCodonFreq[22];
	mQ[1124] = mCodonFreq[26];
	mQ[1132] = mCodonFreq[34];
	mQ[1148] = mCodonFreq[50];
	mQ[1165] = aK*mCodonFreq[6];
	mQ[1174] = aK*mCodonFreq[15];
	mQ[1176] = mCodonFreq[17];
	mQ[1177] = mCodonFreq[18];
	mQ[1178] = -aK*mCodonFreq[6]-aK*mCodonFreq[15]-mCodonFreq[17]-mCodonFreq[18]-aK*mCodonFreq[20]-mCodonFreq[23]-mCodonFreq[27]-mCodonFreq[35]-mCodonFreq[51];
	mQ[1179] = aK*mCodonFreq[20];
	mQ[1182] = mCodonFreq[23];
	mQ[1186] = mCodonFreq[27];
	mQ[1194] = mCodonFreq[35];
	mQ[1210] = mCodonFreq[51];
	mQ[1227] = aK*mCodonFreq[7];
	mQ[1236] = aK*mCodonFreq[16];
	mQ[1237] = mCodonFreq[17];
	mQ[1238] = mCodonFreq[18];
	mQ[1239] = aK*mCodonFreq[19];
	mQ[1240] = -aK*mCodonFreq[7]-aK*mCodonFreq[16]-mCodonFreq[17]-mCodonFreq[18]-aK*mCodonFreq[19]-mCodonFreq[24]-mCodonFreq[28]-mCodonFreq[36]-mCodonFreq[52];
	mQ[1244] = mCodonFreq[24];
	mQ[1248] = mCodonFreq[28];
	mQ[1256] = mCodonFreq[36];
	mQ[1272] = mCodonFreq[52];
	mQ[1289] = aK*mCodonFreq[8];
	mQ[1294] = mCodonFreq[13];
	mQ[1298] = mCodonFreq[17];
	mQ[1302] = -aK*mCodonFreq[8]-mCodonFreq[13]-mCodonFreq[17]-aK*mCodonFreq[22]-mCodonFreq[23]-mCodonFreq[24]-aK*mCodonFreq[25]-mCodonFreq[37]-mCodonFreq[53];
	mQ[1303] = aK*mCodonFreq[22];
	mQ[1304] = mCodonFreq[23];
	mQ[1305] = mCodonFreq[24];
	mQ[1306] = aK*mCodonFreq[25];
	mQ[1318] = mCodonFreq[37];
	mQ[1334] = mCodonFreq[53];
	mQ[1351] = aK*mCodonFreq[9];
	mQ[1356] = mCodonFreq[14];
	mQ[1360] = mCodonFreq[18];
	mQ[1363] = aK*mCodonFreq[21];
	mQ[1364] = -aK*mCodonFreq[9]-mCodonFreq[14]-mCodonFreq[18]-aK*mCodonFreq[21]-mCodonFreq[23]-mCodonFreq[24]-aK*mCodonFreq[26]-mCodonFreq[38]-mCodonFreq[54];
	mQ[1365] = mCodonFreq[23];
	mQ[1366] = mCodonFreq[24];
	mQ[1368] = aK*mCodonFreq[26];
	mQ[1380] = mCodonFreq[38];
	mQ[1396] = mCodonFreq[54];
	mQ[1418] = mCodonFreq[15];
	mQ[1422] = mCodonFreq[19];
	mQ[1424] = mCodonFreq[21];
	mQ[1425] = mCodonFreq[22];
	mQ[1426] = -mCodonFreq[15]-mCodonFreq[19]-mCodonFreq[21]-mCodonFreq[22]-aK*mCodonFreq[24]-aK*mCodonFreq[27]-mCodonFreq[39]-mCodonFreq[55];
	mQ[1427] = aK*mCodonFreq[24];
	mQ[1430] = aK*mCodonFreq[27];
	mQ[1442] = mCodonFreq[39];
	mQ[1458] = mCodonFreq[55];
	mQ[1480] = mCodonFreq[16];
	mQ[1484] = mCodonFreq[20];
	mQ[1485] = mCodonFreq[21];
	mQ[1486] = mCodonFreq[22];
	mQ[1487] = aK*mCodonFreq[23];
	mQ[1488] = -mCodonFreq[16]-mCodonFreq[20]-mCodonFreq[21]-mCodonFreq[22]-aK*mCodonFreq[23]-aK*mCodonFreq[28]-mCodonFreq[40]-mCodonFreq[56];
	mQ[1492] = aK*mCodonFreq[28];
	mQ[1504] = mCodonFreq[40];
	mQ[1520] = mCodonFreq[56];
	mQ[1535] = aK*mCodonFreq[10];
	mQ[1538] = mCodonFreq[13];
	mQ[1542] = mCodonFreq[17];
	mQ[1546] = aK*mCodonFreq[21];
	mQ[1550] = -aK*mCodonFreq[10]-mCodonFreq[13]-mCodonFreq[17]-aK*mCodonFreq[21]-aK*mCodonFreq[26]-mCodonFreq[27]-mCodonFreq[28]-mCodonFreq[41]-mCodonFreq[57];
	mQ[1551] = aK*mCodonFreq[26];
	mQ[1552] = mCodonFreq[27];
	mQ[1553] = mCodonFreq[28];
	mQ[1566] = mCodonFreq[41];
	mQ[1582] = mCodonFreq[57];
	mQ[1597] = aK*mCodonFreq[11];
	mQ[1600] = mCodonFreq[14];
	mQ[1604] = mCodonFreq[18];
	mQ[1608] = aK*mCodonFreq[22];
	mQ[1611] = aK*mCodonFreq[25];
	mQ[1612] = -aK*mCodonFreq[11]-mCodonFreq[14]-mCodonFreq[18]-aK*mCodonFreq[22]-aK*mCodonFreq[25]-mCodonFreq[27]-mCodonFreq[28]-mCodonFreq[42]-mCodonFreq[58];
	mQ[1613] = mCodonFreq[27];
	mQ[1614] = mCodonFreq[28];
	mQ[1628] = mCodonFreq[42];
	mQ[1644] = mCodonFreq[58];
	mQ[1662] = mCodonFreq[15];
	mQ[1666] = mCodonFreq[19];
	mQ[1670] = aK*mCodonFreq[23];
	mQ[1672] = mCodonFreq[25];
	mQ[1673] = mCodonFreq[26];
	mQ[1674] = -mCodonFreq[15]-mCodonFreq[19]-aK*mCodonFreq[23]-mCodonFreq[25]-mCodonFreq[26]-aK*mCodonFreq[28]-mCodonFreq[43]-mCodonFreq[59];
	mQ[1675] = aK*mCodonFreq[28];
	mQ[1690] = mCodonFreq[43];
	mQ[1706] = mCodonFreq[59];
	mQ[1720] = aK*mCodonFreq[12];
	mQ[1724] = mCodonFreq[16];
	mQ[1728] = mCodonFreq[20];
	mQ[1732] = aK*mCodonFreq[24];
	mQ[1733] = mCodonFreq[25];
	mQ[1734] = mCodonFreq[26];
	mQ[1735] = aK*mCodonFreq[27];
	mQ[1736] = -aK*mCodonFreq[12]-mCodonFreq[16]-mCodonFreq[20]-aK*mCodonFreq[24]-mCodonFreq[25]-mCodonFreq[26]-aK*mCodonFreq[27]-mCodonFreq[44]-mCodonFreq[60];
	mQ[1752] = mCodonFreq[44];
	mQ[1768] = mCodonFreq[60];
	mQ[1769] = mCodonFreq[0];
	mQ[1782] = mCodonFreq[13];
	mQ[1798] = -mCodonFreq[0]-mCodonFreq[13]-aK*mCodonFreq[30]-mCodonFreq[31]-mCodonFreq[32]-aK*mCodonFreq[33]-mCodonFreq[37]-mCodonFreq[41]-aK*mCodonFreq[45];
	mQ[1799] = aK*mCodonFreq[30];
	mQ[1800] = mCodonFreq[31];
	mQ[1801] = mCodonFreq[32];
	mQ[1802] = aK*mCodonFreq[33];
	mQ[1806] = mCodonFreq[37];
	mQ[1810] = mCodonFreq[41];
	mQ[1814] = aK*mCodonFreq[45];
	mQ[1831] = mCodonFreq[1];
	mQ[1844] = mCodonFreq[14];
	mQ[1859] = aK*mCodonFreq[29];
	mQ[1860] = -mCodonFreq[1]-mCodonFreq[14]-aK*mCodonFreq[29]-mCodonFreq[31]-mCodonFreq[32]-aK*mCodonFreq[34]-mCodonFreq[38]-mCodonFreq[42]-aK*mCodonFreq[46];
	mQ[1861] = mCodonFreq[31];
	mQ[1862] = mCodonFreq[32];
	mQ[1864] = aK*mCodonFreq[34];
	mQ[1868] = mCodonFreq[38];
	mQ[1872] = mCodonFreq[42];
	mQ[1876] = aK*mCodonFreq[46];
	mQ[1893] = mCodonFreq[2];
	mQ[1906] = mCodonFreq[15];
	mQ[1920] = mCodonFreq[29];
	mQ[1921] = mCodonFreq[30];
	mQ[1922] = -mCodonFreq[2]-mCodonFreq[15]-mCodonFreq[29]-mCodonFreq[30]-aK*mCodonFreq[32]-aK*mCodonFreq[35]-mCodonFreq[39]-mCodonFreq[43]-aK*mCodonFreq[47];
	mQ[1923] = aK*mCodonFreq[32];
	mQ[1926] = aK*mCodonFreq[35];
	mQ[1930] = mCodonFreq[39];
	mQ[1934] = mCodonFreq[43];
	mQ[1938] = aK*mCodonFreq[47];
	mQ[1955] = mCodonFreq[3];
	mQ[1968] = mCodonFreq[16];
	mQ[1981] = mCodonFreq[29];
	mQ[1982] = mCodonFreq[30];
	mQ[1983] = aK*mCodonFreq[31];
	mQ[1984] = -mCodonFreq[3]-mCodonFreq[16]-mCodonFreq[29]-mCodonFreq[30]-aK*mCodonFreq[31]-aK*mCodonFreq[36]-mCodonFreq[40]-mCodonFreq[44]-aK*mCodonFreq[48];
	mQ[1988] = aK*mCodonFreq[36];
	mQ[1992] = mCodonFreq[40];
	mQ[1996] = mCodonFreq[44];
	mQ[2000] = aK*mCodonFreq[48];
	mQ[2017] = mCodonFreq[4];
	mQ[2030] = mCodonFreq[17];
	mQ[2042] = aK*mCodonFreq[29];
	mQ[2046] = -mCodonFreq[4]-mCodonFreq[17]-aK*mCodonFreq[29]-aK*mCodonFreq[34]-mCodonFreq[35]-mCodonFreq[36]-mCodonFreq[37]-mCodonFreq[41]-aK*mCodonFreq[49];
	mQ[2047] = aK*mCodonFreq[34];
	mQ[2048] = mCodonFreq[35];
	mQ[2049] = mCodonFreq[36];
	mQ[2050] = mCodonFreq[37];
	mQ[2054] = mCodonFreq[41];
	mQ[2062] = aK*mCodonFreq[49];
	mQ[2079] = mCodonFreq[5];
	mQ[2092] = mCodonFreq[18];
	mQ[2104] = aK*mCodonFreq[30];
	mQ[2107] = aK*mCodonFreq[33];
	mQ[2108] = -mCodonFreq[5]-mCodonFreq[18]-aK*mCodonFreq[30]-aK*mCodonFreq[33]-mCodonFreq[35]-mCodonFreq[36]-mCodonFreq[38]-mCodonFreq[42]-aK*mCodonFreq[50];
	mQ[2109] = mCodonFreq[35];
	mQ[2110] = mCodonFreq[36];
	mQ[2112] = mCodonFreq[38];
	mQ[2116] = mCodonFreq[42];
	mQ[2124] = aK*mCodonFreq[50];
	mQ[2141] = mCodonFreq[6];
	mQ[2154] = mCodonFreq[19];
	mQ[2166] = aK*mCodonFreq[31];
	mQ[2168] = mCodonFreq[33];
	mQ[2169] = mCodonFreq[34];
	mQ[2170] = -mCodonFreq[6]-mCodonFreq[19]-aK*mCodonFreq[31]-mCodonFreq[33]-mCodonFreq[34]-aK*mCodonFreq[36]-mCodonFreq[39]-mCodonFreq[43]-aK*mCodonFreq[51];
	mQ[2171] = aK*mCodonFreq[36];
	mQ[2174] = mCodonFreq[39];
	mQ[2178] = mCodonFreq[43];
	mQ[2186] = aK*mCodonFreq[51];
	mQ[2203] = mCodonFreq[7];
	mQ[2216] = mCodonFreq[20];
	mQ[2228] = aK*mCodonFreq[32];
	mQ[2229] = mCodonFreq[33];
	mQ[2230] = mCodonFreq[34];
	mQ[2231] = aK*mCodonFreq[35];
	mQ[2232] = -mCodonFreq[7]-mCodonFreq[20]-aK*mCodonFreq[32]-mCodonFreq[33]-mCodonFreq[34]-aK*mCodonFreq[35]-mCodonFreq[40]-mCodonFreq[44]-aK*mCodonFreq[52];
	mQ[2236] = mCodonFreq[40];
	mQ[2240] = mCodonFreq[44];
	mQ[2248] = aK*mCodonFreq[52];
	mQ[2265] = mCodonFreq[8];
	mQ[2278] = mCodonFreq[21];
	mQ[2286] = mCodonFreq[29];
	mQ[2290] = mCodonFreq[33];
	mQ[2294] = -mCodonFreq[8]-mCodonFreq[21]-mCodonFreq[29]-mCodonFreq[33]-aK*mCodonFreq[38]-mCodonFreq[39]-mCodonFreq[40]-aK*mCodonFreq[41]-aK*mCodonFreq[53];
	mQ[2295] = aK*mCodonFreq[38];
	mQ[2296] = mCodonFreq[39];
	mQ[2297] = mCodonFreq[40];
	mQ[2298] = aK*mCodonFreq[41];
	mQ[2310] = aK*mCodonFreq[53];
	mQ[2327] = mCodonFreq[9];
	mQ[2340] = mCodonFreq[22];
	mQ[2348] = mCodonFreq[30];
	mQ[2352] = mCodonFreq[34];
	mQ[2355] = aK*mCodonFreq[37];
	mQ[2356] = -mCodonFreq[9]-mCodonFreq[22]-mCodonFreq[30]-mCodonFreq[34]-aK*mCodonFreq[37]-mCodonFreq[39]-mCodonFreq[40]-aK*mCodonFreq[42]-aK*mCodonFreq[54];
	mQ[2357] = mCodonFreq[39];
	mQ[2358] = mCodonFreq[40];
	mQ[2360] = aK*mCodonFreq[42];
	mQ[2372] = aK*mCodonFreq[54];
	mQ[2402] = mCodonFreq[23];
	mQ[2410] = mCodonFreq[31];
	mQ[2414] = mCodonFreq[35];
	mQ[2416] = mCodonFreq[37];
	mQ[2417] = mCodonFreq[38];
	mQ[2418] = -mCodonFreq[23]-mCodonFreq[31]-mCodonFreq[35]-mCodonFreq[37]-mCodonFreq[38]-aK*mCodonFreq[40]-aK*mCodonFreq[43]-aK*mCodonFreq[55];
	mQ[2419] = aK*mCodonFreq[40];
	mQ[2422] = aK*mCodonFreq[43];
	mQ[2434] = aK*mCodonFreq[55];
	mQ[2464] = mCodonFreq[24];
	mQ[2472] = mCodonFreq[32];
	mQ[2476] = mCodonFreq[36];
	mQ[2477] = mCodonFreq[37];
	mQ[2478] = mCodonFreq[38];
	mQ[2479] = aK*mCodonFreq[39];
	mQ[2480] = -mCodonFreq[24]-mCodonFreq[32]-mCodonFreq[36]-mCodonFreq[37]-mCodonFreq[38]-aK*mCodonFreq[39]-aK*mCodonFreq[44]-aK*mCodonFreq[56];
	mQ[2484] = aK*mCodonFreq[44];
	mQ[2496] = aK*mCodonFreq[56];
	mQ[2511] = mCodonFreq[10];
	mQ[2526] = mCodonFreq[25];
	mQ[2530] = mCodonFreq[29];
	mQ[2534] = mCodonFreq[33];
	mQ[2538] = aK*mCodonFreq[37];
	mQ[2542] = -mCodonFreq[10]-mCodonFreq[25]-mCodonFreq[29]-mCodonFreq[33]-aK*mCodonFreq[37]-aK*mCodonFreq[42]-mCodonFreq[43]-mCodonFreq[44]-aK*mCodonFreq[57];
	mQ[2543] = aK*mCodonFreq[42];
	mQ[2544] = mCodonFreq[43];
	mQ[2545] = mCodonFreq[44];
	mQ[2558] = aK*mCodonFreq[57];
	mQ[2573] = mCodonFreq[11];
	mQ[2588] = mCodonFreq[26];
	mQ[2592] = mCodonFreq[30];
	mQ[2596] = mCodonFreq[34];
	mQ[2600] = aK*mCodonFreq[38];
	mQ[2603] = aK*mCodonFreq[41];
	mQ[2604] = -mCodonFreq[11]-mCodonFreq[26]-mCodonFreq[30]-mCodonFreq[34]-aK*mCodonFreq[38]-aK*mCodonFreq[41]-mCodonFreq[43]-mCodonFreq[44]-aK*mCodonFreq[58];
	mQ[2605] = mCodonFreq[43];
	mQ[2606] = mCodonFreq[44];
	mQ[2620] = aK*mCodonFreq[58];
	mQ[2650] = mCodonFreq[27];
	mQ[2654] = mCodonFreq[31];
	mQ[2658] = mCodonFreq[35];
	mQ[2662] = aK*mCodonFreq[39];
	mQ[2664] = mCodonFreq[41];
	mQ[2665] = mCodonFreq[42];
	mQ[2666] = -mCodonFreq[27]-mCodonFreq[31]-mCodonFreq[35]-aK*mCodonFreq[39]-mCodonFreq[41]-mCodonFreq[42]-aK*mCodonFreq[44]-aK*mCodonFreq[59];
	mQ[2667] = aK*mCodonFreq[44];
	mQ[2682] = aK*mCodonFreq[59];
	mQ[2696] = mCodonFreq[12];
	mQ[2712] = mCodonFreq[28];
	mQ[2716] = mCodonFreq[32];
	mQ[2720] = mCodonFreq[36];
	mQ[2724] = aK*mCodonFreq[40];
	mQ[2725] = mCodonFreq[41];
	mQ[2726] = mCodonFreq[42];
	mQ[2727] = aK*mCodonFreq[43];
	mQ[2728] = -mCodonFreq[12]-mCodonFreq[28]-mCodonFreq[32]-mCodonFreq[36]-aK*mCodonFreq[40]-mCodonFreq[41]-mCodonFreq[42]-aK*mCodonFreq[43]-aK*mCodonFreq[60];
	mQ[2744] = aK*mCodonFreq[60];
	mQ[2745] = mCodonFreq[0];
	mQ[2758] = mCodonFreq[13];
	mQ[2774] = aK*mCodonFreq[29];
	mQ[2790] = -mCodonFreq[0]-mCodonFreq[13]-aK*mCodonFreq[29]-aK*mCodonFreq[46]-mCodonFreq[47]-mCodonFreq[48]-aK*mCodonFreq[49]-mCodonFreq[53]-mCodonFreq[57];
	mQ[2791] = aK*mCodonFreq[46];
	mQ[2792] = mCodonFreq[47];
	mQ[2793] = mCodonFreq[48];
	mQ[2794] = aK*mCodonFreq[49];
	mQ[2798] = mCodonFreq[53];
	mQ[2802] = mCodonFreq[57];
	mQ[2807] = mCodonFreq[1];
	mQ[2820] = mCodonFreq[14];
	mQ[2836] = aK*mCodonFreq[30];
	mQ[2851] = aK*mCodonFreq[45];
	mQ[2852] = -mCodonFreq[1]-mCodonFreq[14]-aK*mCodonFreq[30]-aK*mCodonFreq[45]-mCodonFreq[47]-mCodonFreq[48]-aK*mCodonFreq[50]-mCodonFreq[54]-mCodonFreq[58];
	mQ[2853] = mCodonFreq[47];
	mQ[2854] = mCodonFreq[48];
	mQ[2856] = aK*mCodonFreq[50];
	mQ[2860] = mCodonFreq[54];
	mQ[2864] = mCodonFreq[58];
	mQ[2869] = mCodonFreq[2];
	mQ[2882] = mCodonFreq[15];
	mQ[2898] = aK*mCodonFreq[31];
	mQ[2912] = mCodonFreq[45];
	mQ[2913] = mCodonFreq[46];
	mQ[2914] = -mCodonFreq[2]-mCodonFreq[15]-aK*mCodonFreq[31]-mCodonFreq[45]-mCodonFreq[46]-aK*mCodonFreq[48]-aK*mCodonFreq[51]-mCodonFreq[55]-mCodonFreq[59];
	mQ[2915] = aK*mCodonFreq[48];
	mQ[2918] = aK*mCodonFreq[51];
	mQ[2922] = mCodonFreq[55];
	mQ[2926] = mCodonFreq[59];
	mQ[2931] = mCodonFreq[3];
	mQ[2944] = mCodonFreq[16];
	mQ[2960] = aK*mCodonFreq[32];
	mQ[2973] = mCodonFreq[45];
	mQ[2974] = mCodonFreq[46];
	mQ[2975] = aK*mCodonFreq[47];
	mQ[2976] = -mCodonFreq[3]-mCodonFreq[16]-aK*mCodonFreq[32]-mCodonFreq[45]-mCodonFreq[46]-aK*mCodonFreq[47]-aK*mCodonFreq[52]-mCodonFreq[56]-mCodonFreq[60];
	mQ[2980] = aK*mCodonFreq[52];
	mQ[2984] = mCodonFreq[56];
	mQ[2988] = mCodonFreq[60];
	mQ[2993] = mCodonFreq[4];
	mQ[3006] = mCodonFreq[17];
	mQ[3022] = aK*mCodonFreq[33];
	mQ[3034] = aK*mCodonFreq[45];
	mQ[3038] = -mCodonFreq[4]-mCodonFreq[17]-aK*mCodonFreq[33]-aK*mCodonFreq[45]-aK*mCodonFreq[50]-mCodonFreq[51]-mCodonFreq[52]-mCodonFreq[53]-mCodonFreq[57];
	mQ[3039] = aK*mCodonFreq[50];
	mQ[3040] = mCodonFreq[51];
	mQ[3041] = mCodonFreq[52];
	mQ[3042] = mCodonFreq[53];
	mQ[3046] = mCodonFreq[57];
	mQ[3055] = mCodonFreq[5];
	mQ[3068] = mCodonFreq[18];
	mQ[3084] = aK*mCodonFreq[34];
	mQ[3096] = aK*mCodonFreq[46];
	mQ[3099] = aK*mCodonFreq[49];
	mQ[3100] = -mCodonFreq[5]-mCodonFreq[18]-aK*mCodonFreq[34]-aK*mCodonFreq[46]-aK*mCodonFreq[49]-mCodonFreq[51]-mCodonFreq[52]-mCodonFreq[54]-mCodonFreq[58];
	mQ[3101] = mCodonFreq[51];
	mQ[3102] = mCodonFreq[52];
	mQ[3104] = mCodonFreq[54];
	mQ[3108] = mCodonFreq[58];
	mQ[3117] = mCodonFreq[6];
	mQ[3130] = mCodonFreq[19];
	mQ[3146] = aK*mCodonFreq[35];
	mQ[3158] = aK*mCodonFreq[47];
	mQ[3160] = mCodonFreq[49];
	mQ[3161] = mCodonFreq[50];
	mQ[3162] = -mCodonFreq[6]-mCodonFreq[19]-aK*mCodonFreq[35]-aK*mCodonFreq[47]-mCodonFreq[49]-mCodonFreq[50]-aK*mCodonFreq[52]-mCodonFreq[55]-mCodonFreq[59];
	mQ[3163] = aK*mCodonFreq[52];
	mQ[3166] = mCodonFreq[55];
	mQ[3170] = mCodonFreq[59];
	mQ[3179] = mCodonFreq[7];
	mQ[3192] = mCodonFreq[20];
	mQ[3208] = aK*mCodonFreq[36];
	mQ[3220] = aK*mCodonFreq[48];
	mQ[3221] = mCodonFreq[49];
	mQ[3222] = mCodonFreq[50];
	mQ[3223] = aK*mCodonFreq[51];
	mQ[3224] = -mCodonFreq[7]-mCodonFreq[20]-aK*mCodonFreq[36]-aK*mCodonFreq[48]-mCodonFreq[49]-mCodonFreq[50]-aK*mCodonFreq[51]-mCodonFreq[56]-mCodonFreq[60];
	mQ[3228] = mCodonFreq[56];
	mQ[3232] = mCodonFreq[60];
	mQ[3241] = mCodonFreq[8];
	mQ[3254] = mCodonFreq[21];
	mQ[3270] = aK*mCodonFreq[37];
	mQ[3278] = mCodonFreq[45];
	mQ[3282] = mCodonFreq[49];
	mQ[3286] = -mCodonFreq[8]-mCodonFreq[21]-aK*mCodonFreq[37]-mCodonFreq[45]-mCodonFreq[49]-aK*mCodonFreq[54]-mCodonFreq[55]-mCodonFreq[56]-aK*mCodonFreq[57];
	mQ[3287] = aK*mCodonFreq[54];
	mQ[3288] = mCodonFreq[55];
	mQ[3289] = mCodonFreq[56];
	mQ[3290] = aK*mCodonFreq[57];
	mQ[3303] = mCodonFreq[9];
	mQ[3316] = mCodonFreq[22];
	mQ[3332] = aK*mCodonFreq[38];
	mQ[3340] = mCodonFreq[46];
	mQ[3344] = mCodonFreq[50];
	mQ[3347] = aK*mCodonFreq[53];
	mQ[3348] = -mCodonFreq[9]-mCodonFreq[22]-aK*mCodonFreq[38]-mCodonFreq[46]-mCodonFreq[50]-aK*mCodonFreq[53]-mCodonFreq[55]-mCodonFreq[56]-aK*mCodonFreq[58];
	mQ[3349] = mCodonFreq[55];
	mQ[3350] = mCodonFreq[56];
	mQ[3352] = aK*mCodonFreq[58];
	mQ[3378] = mCodonFreq[23];
	mQ[3394] = aK*mCodonFreq[39];
	mQ[3402] = mCodonFreq[47];
	mQ[3406] = mCodonFreq[51];
	mQ[3408] = mCodonFreq[53];
	mQ[3409] = mCodonFreq[54];
	mQ[3410] = -mCodonFreq[23]-aK*mCodonFreq[39]-mCodonFreq[47]-mCodonFreq[51]-mCodonFreq[53]-mCodonFreq[54]-aK*mCodonFreq[56]-aK*mCodonFreq[59];
	mQ[3411] = aK*mCodonFreq[56];
	mQ[3414] = aK*mCodonFreq[59];
	mQ[3440] = mCodonFreq[24];
	mQ[3456] = aK*mCodonFreq[40];
	mQ[3464] = mCodonFreq[48];
	mQ[3468] = mCodonFreq[52];
	mQ[3469] = mCodonFreq[53];
	mQ[3470] = mCodonFreq[54];
	mQ[3471] = aK*mCodonFreq[55];
	mQ[3472] = -mCodonFreq[24]-aK*mCodonFreq[40]-mCodonFreq[48]-mCodonFreq[52]-mCodonFreq[53]-mCodonFreq[54]-aK*mCodonFreq[55]-aK*mCodonFreq[60];
	mQ[3476] = aK*mCodonFreq[60];
	mQ[3487] = mCodonFreq[10];
	mQ[3502] = mCodonFreq[25];
	mQ[3518] = aK*mCodonFreq[41];
	mQ[3522] = mCodonFreq[45];
	mQ[3526] = mCodonFreq[49];
	mQ[3530] = aK*mCodonFreq[53];
	mQ[3534] = -mCodonFreq[10]-mCodonFreq[25]-aK*mCodonFreq[41]-mCodonFreq[45]-mCodonFreq[49]-aK*mCodonFreq[53]-aK*mCodonFreq[58]-mCodonFreq[59]-mCodonFreq[60];
	mQ[3535] = aK*mCodonFreq[58];
	mQ[3536] = mCodonFreq[59];
	mQ[3537] = mCodonFreq[60];
	mQ[3549] = mCodonFreq[11];
	mQ[3564] = mCodonFreq[26];
	mQ[3580] = aK*mCodonFreq[42];
	mQ[3584] = mCodonFreq[46];
	mQ[3588] = mCodonFreq[50];
	mQ[3592] = aK*mCodonFreq[54];
	mQ[3595] = aK*mCodonFreq[57];
	mQ[3596] = -mCodonFreq[11]-mCodonFreq[26]-aK*mCodonFreq[42]-mCodonFreq[46]-mCodonFreq[50]-aK*mCodonFreq[54]-aK*mCodonFreq[57]-mCodonFreq[59]-mCodonFreq[60];
	mQ[3597] = mCodonFreq[59];
	mQ[3598] = mCodonFreq[60];
	mQ[3626] = mCodonFreq[27];
	mQ[3642] = aK*mCodonFreq[43];
	mQ[3646] = mCodonFreq[47];
	mQ[3650] = mCodonFreq[51];
	mQ[3654] = aK*mCodonFreq[55];
	mQ[3656] = mCodonFreq[57];
	mQ[3657] = mCodonFreq[58];
	mQ[3658] = -mCodonFreq[27]-aK*mCodonFreq[43]-mCodonFreq[47]-mCodonFreq[51]-aK*mCodonFreq[55]-mCodonFreq[57]-mCodonFreq[58]-aK*mCodonFreq[60];
	mQ[3659] = aK*mCodonFreq[60];
	mQ[3672] = mCodonFreq[12];
	mQ[3688] = mCodonFreq[28];
	mQ[3704] = aK*mCodonFreq[44];
	mQ[3708] = mCodonFreq[48];
	mQ[3712] = mCodonFreq[52];
	mQ[3716] = aK*mCodonFreq[56];
	mQ[3717] = mCodonFreq[57];
	mQ[3718] = mCodonFreq[58];
	mQ[3719] = aK*mCodonFreq[59];
	mQ[3720] = -mCodonFreq[12]-mCodonFreq[28]-aK*mCodonFreq[44]-mCodonFreq[48]-mCodonFreq[52]-aK*mCodonFreq[56]-mCodonFreq[57]-mCodonFreq[58]-aK*mCodonFreq[59];

	// Compute the scale factor
	double scale_q =     aK*mCodonFreq[0]*mCodonFreq[1]
	             +   mCodonFreq[0]*mCodonFreq[2]
	             +   mCodonFreq[0]*mCodonFreq[3]
	             + aK*mCodonFreq[0]*mCodonFreq[4]
	             +   mCodonFreq[0]*mCodonFreq[8]
	             +   mCodonFreq[0]*mCodonFreq[10]
	             + aK*mCodonFreq[0]*mCodonFreq[13]
	             +   mCodonFreq[0]*mCodonFreq[29]
	             +   mCodonFreq[0]*mCodonFreq[45]
	             +   mCodonFreq[1]*mCodonFreq[2]
	             +   mCodonFreq[1]*mCodonFreq[3]
	             + aK*mCodonFreq[1]*mCodonFreq[5]
	             +   mCodonFreq[1]*mCodonFreq[9]
	             +   mCodonFreq[1]*mCodonFreq[11]
	             + aK*mCodonFreq[1]*mCodonFreq[14]
	             +   mCodonFreq[1]*mCodonFreq[30]
	             +   mCodonFreq[1]*mCodonFreq[46]
	             + aK*mCodonFreq[2]*mCodonFreq[3]
	             + aK*mCodonFreq[2]*mCodonFreq[6]
	             + aK*mCodonFreq[2]*mCodonFreq[15]
	             +   mCodonFreq[2]*mCodonFreq[31]
	             +   mCodonFreq[2]*mCodonFreq[47]
	             + aK*mCodonFreq[3]*mCodonFreq[7]
	             +   mCodonFreq[3]*mCodonFreq[12]
	             + aK*mCodonFreq[3]*mCodonFreq[16]
	             +   mCodonFreq[3]*mCodonFreq[32]
	             +   mCodonFreq[3]*mCodonFreq[48]
	             + aK*mCodonFreq[4]*mCodonFreq[5]
	             +   mCodonFreq[4]*mCodonFreq[6]
	             +   mCodonFreq[4]*mCodonFreq[7]
	             +   mCodonFreq[4]*mCodonFreq[8]
	             +   mCodonFreq[4]*mCodonFreq[10]
	             + aK*mCodonFreq[4]*mCodonFreq[17]
	             +   mCodonFreq[4]*mCodonFreq[33]
	             +   mCodonFreq[4]*mCodonFreq[49]
	             +   mCodonFreq[5]*mCodonFreq[6]
	             +   mCodonFreq[5]*mCodonFreq[7]
	             +   mCodonFreq[5]*mCodonFreq[9]
	             +   mCodonFreq[5]*mCodonFreq[11]
	             + aK*mCodonFreq[5]*mCodonFreq[18]
	             +   mCodonFreq[5]*mCodonFreq[34]
	             +   mCodonFreq[5]*mCodonFreq[50]
	             + aK*mCodonFreq[6]*mCodonFreq[7]
	             + aK*mCodonFreq[6]*mCodonFreq[19]
	             +   mCodonFreq[6]*mCodonFreq[35]
	             +   mCodonFreq[6]*mCodonFreq[51]
	             +   mCodonFreq[7]*mCodonFreq[12]
	             + aK*mCodonFreq[7]*mCodonFreq[20]
	             +   mCodonFreq[7]*mCodonFreq[36]
	             +   mCodonFreq[7]*mCodonFreq[52]
	             + aK*mCodonFreq[8]*mCodonFreq[9]
	             + aK*mCodonFreq[8]*mCodonFreq[10]
	             + aK*mCodonFreq[8]*mCodonFreq[21]
	             +   mCodonFreq[8]*mCodonFreq[37]
	             +   mCodonFreq[8]*mCodonFreq[53]
	             + aK*mCodonFreq[9]*mCodonFreq[11]
	             + aK*mCodonFreq[9]*mCodonFreq[22]
	             +   mCodonFreq[9]*mCodonFreq[38]
	             +   mCodonFreq[9]*mCodonFreq[54]
	             + aK*mCodonFreq[10]*mCodonFreq[11]
	             +   mCodonFreq[10]*mCodonFreq[12]
	             + aK*mCodonFreq[10]*mCodonFreq[25]
	             +   mCodonFreq[10]*mCodonFreq[41]
	             +   mCodonFreq[10]*mCodonFreq[57]
	             +   mCodonFreq[11]*mCodonFreq[12]
	             + aK*mCodonFreq[11]*mCodonFreq[26]
	             +   mCodonFreq[11]*mCodonFreq[42]
	             +   mCodonFreq[11]*mCodonFreq[58]
	             + aK*mCodonFreq[12]*mCodonFreq[28]
	             +   mCodonFreq[12]*mCodonFreq[44]
	             +   mCodonFreq[12]*mCodonFreq[60]
	             + aK*mCodonFreq[13]*mCodonFreq[14]
	             +   mCodonFreq[13]*mCodonFreq[15]
	             +   mCodonFreq[13]*mCodonFreq[16]
	             + aK*mCodonFreq[13]*mCodonFreq[17]
	             +   mCodonFreq[13]*mCodonFreq[21]
	             +   mCodonFreq[13]*mCodonFreq[25]
	             +   mCodonFreq[13]*mCodonFreq[29]
	             +   mCodonFreq[13]*mCodonFreq[45]
	             +   mCodonFreq[14]*mCodonFreq[15]
	             +   mCodonFreq[14]*mCodonFreq[16]
	             + aK*mCodonFreq[14]*mCodonFreq[18]
	             +   mCodonFreq[14]*mCodonFreq[22]
	             +   mCodonFreq[14]*mCodonFreq[26]
	             +   mCodonFreq[14]*mCodonFreq[30]
	             +   mCodonFreq[14]*mCodonFreq[46]
	             + aK*mCodonFreq[15]*mCodonFreq[16]
	             + aK*mCodonFreq[15]*mCodonFreq[19]
	             +   mCodonFreq[15]*mCodonFreq[23]
	             +   mCodonFreq[15]*mCodonFreq[27]
	             +   mCodonFreq[15]*mCodonFreq[31]
	             +   mCodonFreq[15]*mCodonFreq[47]
	             + aK*mCodonFreq[16]*mCodonFreq[20]
	             +   mCodonFreq[16]*mCodonFreq[24]
	             +   mCodonFreq[16]*mCodonFreq[28]
	             +   mCodonFreq[16]*mCodonFreq[32]
	             +   mCodonFreq[16]*mCodonFreq[48]
	             + aK*mCodonFreq[17]*mCodonFreq[18]
	             +   mCodonFreq[17]*mCodonFreq[19]
	             +   mCodonFreq[17]*mCodonFreq[20]
	             +   mCodonFreq[17]*mCodonFreq[21]
	             +   mCodonFreq[17]*mCodonFreq[25]
	             +   mCodonFreq[17]*mCodonFreq[33]
	             +   mCodonFreq[17]*mCodonFreq[49]
	             +   mCodonFreq[18]*mCodonFreq[19]
	             +   mCodonFreq[18]*mCodonFreq[20]
	             +   mCodonFreq[18]*mCodonFreq[22]
	             +   mCodonFreq[18]*mCodonFreq[26]
	             +   mCodonFreq[18]*mCodonFreq[34]
	             +   mCodonFreq[18]*mCodonFreq[50]
	             + aK*mCodonFreq[19]*mCodonFreq[20]
	             +   mCodonFreq[19]*mCodonFreq[23]
	             +   mCodonFreq[19]*mCodonFreq[27]
	             +   mCodonFreq[19]*mCodonFreq[35]
	             +   mCodonFreq[19]*mCodonFreq[51]
	             +   mCodonFreq[20]*mCodonFreq[24]
	             +   mCodonFreq[20]*mCodonFreq[28]
	             +   mCodonFreq[20]*mCodonFreq[36]
	             +   mCodonFreq[20]*mCodonFreq[52]
	             + aK*mCodonFreq[21]*mCodonFreq[22]
	             +   mCodonFreq[21]*mCodonFreq[23]
	             +   mCodonFreq[21]*mCodonFreq[24]
	             + aK*mCodonFreq[21]*mCodonFreq[25]
	             +   mCodonFreq[21]*mCodonFreq[37]
	             +   mCodonFreq[21]*mCodonFreq[53]
	             +   mCodonFreq[22]*mCodonFreq[23]
	             +   mCodonFreq[22]*mCodonFreq[24]
	             + aK*mCodonFreq[22]*mCodonFreq[26]
	             +   mCodonFreq[22]*mCodonFreq[38]
	             +   mCodonFreq[22]*mCodonFreq[54]
	             + aK*mCodonFreq[23]*mCodonFreq[24]
	             + aK*mCodonFreq[23]*mCodonFreq[27]
	             +   mCodonFreq[23]*mCodonFreq[39]
	             +   mCodonFreq[23]*mCodonFreq[55]
	             + aK*mCodonFreq[24]*mCodonFreq[28]
	             +   mCodonFreq[24]*mCodonFreq[40]
	             +   mCodonFreq[24]*mCodonFreq[56]
	             + aK*mCodonFreq[25]*mCodonFreq[26]
	             +   mCodonFreq[25]*mCodonFreq[27]
	             +   mCodonFreq[25]*mCodonFreq[28]
	             +   mCodonFreq[25]*mCodonFreq[41]
	             +   mCodonFreq[25]*mCodonFreq[57]
	             +   mCodonFreq[26]*mCodonFreq[27]
	             +   mCodonFreq[26]*mCodonFreq[28]
	             +   mCodonFreq[26]*mCodonFreq[42]
	             +   mCodonFreq[26]*mCodonFreq[58]
	             + aK*mCodonFreq[27]*mCodonFreq[28]
	             +   mCodonFreq[27]*mCodonFreq[43]
	             +   mCodonFreq[27]*mCodonFreq[59]
	             +   mCodonFreq[28]*mCodonFreq[44]
	             +   mCodonFreq[28]*mCodonFreq[60]
	             + aK*mCodonFreq[29]*mCodonFreq[30]
	             +   mCodonFreq[29]*mCodonFreq[31]
	             +   mCodonFreq[29]*mCodonFreq[32]
	             + aK*mCodonFreq[29]*mCodonFreq[33]
	             +   mCodonFreq[29]*mCodonFreq[37]
	             +   mCodonFreq[29]*mCodonFreq[41]
	             + aK*mCodonFreq[29]*mCodonFreq[45]
	             +   mCodonFreq[30]*mCodonFreq[31]
	             +   mCodonFreq[30]*mCodonFreq[32]
	             + aK*mCodonFreq[30]*mCodonFreq[34]
	             +   mCodonFreq[30]*mCodonFreq[38]
	             +   mCodonFreq[30]*mCodonFreq[42]
	             + aK*mCodonFreq[30]*mCodonFreq[46]
	             + aK*mCodonFreq[31]*mCodonFreq[32]
	             + aK*mCodonFreq[31]*mCodonFreq[35]
	             +   mCodonFreq[31]*mCodonFreq[39]
	             +   mCodonFreq[31]*mCodonFreq[43]
	             + aK*mCodonFreq[31]*mCodonFreq[47]
	             + aK*mCodonFreq[32]*mCodonFreq[36]
	             +   mCodonFreq[32]*mCodonFreq[40]
	             +   mCodonFreq[32]*mCodonFreq[44]
	             + aK*mCodonFreq[32]*mCodonFreq[48]
	             + aK*mCodonFreq[33]*mCodonFreq[34]
	             +   mCodonFreq[33]*mCodonFreq[35]
	             +   mCodonFreq[33]*mCodonFreq[36]
	             +   mCodonFreq[33]*mCodonFreq[37]
	             +   mCodonFreq[33]*mCodonFreq[41]
	             + aK*mCodonFreq[33]*mCodonFreq[49]
	             +   mCodonFreq[34]*mCodonFreq[35]
	             +   mCodonFreq[34]*mCodonFreq[36]
	             +   mCodonFreq[34]*mCodonFreq[38]
	             +   mCodonFreq[34]*mCodonFreq[42]
	             + aK*mCodonFreq[34]*mCodonFreq[50]
	             + aK*mCodonFreq[35]*mCodonFreq[36]
	             +   mCodonFreq[35]*mCodonFreq[39]
	             +   mCodonFreq[35]*mCodonFreq[43]
	             + aK*mCodonFreq[35]*mCodonFreq[51]
	             +   mCodonFreq[36]*mCodonFreq[40]
	             +   mCodonFreq[36]*mCodonFreq[44]
	             + aK*mCodonFreq[36]*mCodonFreq[52]
	             + aK*mCodonFreq[37]*mCodonFreq[38]
	             +   mCodonFreq[37]*mCodonFreq[39]
	             +   mCodonFreq[37]*mCodonFreq[40]
	             + aK*mCodonFreq[37]*mCodonFreq[41]
	             + aK*mCodonFreq[37]*mCodonFreq[53]
	             +   mCodonFreq[38]*mCodonFreq[39]
	             +   mCodonFreq[38]*mCodonFreq[40]
	             + aK*mCodonFreq[38]*mCodonFreq[42]
	             + aK*mCodonFreq[38]*mCodonFreq[54]
	             + aK*mCodonFreq[39]*mCodonFreq[40]
	             + aK*mCodonFreq[39]*mCodonFreq[43]
	             + aK*mCodonFreq[39]*mCodonFreq[55]
	             + aK*mCodonFreq[40]*mCodonFreq[44]
	             + aK*mCodonFreq[40]*mCodonFreq[56]
	             + aK*mCodonFreq[41]*mCodonFreq[42]
	             +   mCodonFreq[41]*mCodonFreq[43]
	             +   mCodonFreq[41]*mCodonFreq[44]
	             + aK*mCodonFreq[41]*mCodonFreq[57]
	             +   mCodonFreq[42]*mCodonFreq[43]
	             +   mCodonFreq[42]*mCodonFreq[44]
	             + aK*mCodonFreq[42]*mCodonFreq[58]
	             + aK*mCodonFreq[43]*mCodonFreq[44]
	             + aK*mCodonFreq[43]*mCodonFreq[59]
	             + aK*mCodonFreq[44]*mCodonFreq[60]
	             + aK*mCodonFreq[45]*mCodonFreq[46]
	             +   mCodonFreq[45]*mCodonFreq[47]
	             +   mCodonFreq[45]*mCodonFreq[48]
	             + aK*mCodonFreq[45]*mCodonFreq[49]
	             +   mCodonFreq[45]*mCodonFreq[53]
	             +   mCodonFreq[45]*mCodonFreq[57]
	             +   mCodonFreq[46]*mCodonFreq[47]
	             +   mCodonFreq[46]*mCodonFreq[48]
	             + aK*mCodonFreq[46]*mCodonFreq[50]
	             +   mCodonFreq[46]*mCodonFreq[54]
	             +   mCodonFreq[46]*mCodonFreq[58]
	             + aK*mCodonFreq[47]*mCodonFreq[48]
	             + aK*mCodonFreq[47]*mCodonFreq[51]
	             +   mCodonFreq[47]*mCodonFreq[55]
	             +   mCodonFreq[47]*mCodonFreq[59]
	             + aK*mCodonFreq[48]*mCodonFreq[52]
	             +   mCodonFreq[48]*mCodonFreq[56]
	             +   mCodonFreq[48]*mCodonFreq[60]
	             + aK*mCodonFreq[49]*mCodonFreq[50]
	             +   mCodonFreq[49]*mCodonFreq[51]
	             +   mCodonFreq[49]*mCodonFreq[52]
	             +   mCodonFreq[49]*mCodonFreq[53]
	             +   mCodonFreq[49]*mCodonFreq[57]
	             +   mCodonFreq[50]*mCodonFreq[51]
	             +   mCodonFreq[50]*mCodonFreq[52]
	             +   mCodonFreq[50]*mCodonFreq[54]
	             +   mCodonFreq[50]*mCodonFreq[58]
	             + aK*mCodonFreq[51]*mCodonFreq[52]
	             +   mCodonFreq[51]*mCodonFreq[55]
	             +   mCodonFreq[51]*mCodonFreq[59]
	             +   mCodonFreq[52]*mCodonFreq[56]
	             +   mCodonFreq[52]*mCodonFreq[60]
	             + aK*mCodonFreq[53]*mCodonFreq[54]
	             +   mCodonFreq[53]*mCodonFreq[55]
	             +   mCodonFreq[53]*mCodonFreq[56]
	             + aK*mCodonFreq[53]*mCodonFreq[57]
	             +   mCodonFreq[54]*mCodonFreq[55]
	             +   mCodonFreq[54]*mCodonFreq[56]
	             + aK*mCodonFreq[54]*mCodonFreq[58]
	             + aK*mCodonFreq[55]*mCodonFreq[56]
	             + aK*mCodonFreq[55]*mCodonFreq[59]
	             + aK*mCodonFreq[56]*mCodonFreq[60]
	             + aK*mCodonFreq[57]*mCodonFreq[58]
	             +   mCodonFreq[57]*mCodonFreq[59]
	             +   mCodonFreq[57]*mCodonFreq[60]
	             +   mCodonFreq[58]*mCodonFreq[59]
	             +   mCodonFreq[58]*mCodonFreq[60]
	             + aK*mCodonFreq[59]*mCodonFreq[60];

	return scale_q*2.;
}

