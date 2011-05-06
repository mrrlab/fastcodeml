// Automatically generated file using genQ.r script, don't modify!

#include "TransitionMatrix.h"

double TransitionMatrix::fillQ(double aOmega, double aK, const double* aCodonFreq)
{
	mQ[   0] = -aK*aCodonFreq[1]-aOmega*aCodonFreq[2]-aOmega*aCodonFreq[3]-aK*aOmega*aCodonFreq[4]-aOmega*aCodonFreq[8]-aOmega*aCodonFreq[10]-aK*aOmega*aCodonFreq[13]-aOmega*aCodonFreq[29]-aOmega*aCodonFreq[45];
	mQ[   1] = aK*aCodonFreq[1];
	mQ[   2] = aOmega*aCodonFreq[2];
	mQ[   3] = aOmega*aCodonFreq[3];
	mQ[   4] = aK*aOmega*aCodonFreq[4];
	mQ[   8] = aOmega*aCodonFreq[8];
	mQ[  10] = aOmega*aCodonFreq[10];
	mQ[  13] = aK*aOmega*aCodonFreq[13];
	mQ[  29] = aOmega*aCodonFreq[29];
	mQ[  45] = aOmega*aCodonFreq[45];
	mQ[  61] = aK*aCodonFreq[0];
	mQ[  62] = -aK*aCodonFreq[0]-aOmega*aCodonFreq[2]-aOmega*aCodonFreq[3]-aK*aOmega*aCodonFreq[5]-aOmega*aCodonFreq[9]-aOmega*aCodonFreq[11]-aK*aOmega*aCodonFreq[14]-aOmega*aCodonFreq[30]-aOmega*aCodonFreq[46];
	mQ[  63] = aOmega*aCodonFreq[2];
	mQ[  64] = aOmega*aCodonFreq[3];
	mQ[  66] = aK*aOmega*aCodonFreq[5];
	mQ[  70] = aOmega*aCodonFreq[9];
	mQ[  72] = aOmega*aCodonFreq[11];
	mQ[  75] = aK*aOmega*aCodonFreq[14];
	mQ[  91] = aOmega*aCodonFreq[30];
	mQ[ 107] = aOmega*aCodonFreq[46];
	mQ[ 122] = aOmega*aCodonFreq[0];
	mQ[ 123] = aOmega*aCodonFreq[1];
	mQ[ 124] = -aOmega*aCodonFreq[0]-aOmega*aCodonFreq[1]-aK*aCodonFreq[3]-aK*aOmega*aCodonFreq[6]-aK*aCodonFreq[15]-aOmega*aCodonFreq[31]-aOmega*aCodonFreq[47];
	mQ[ 125] = aK*aCodonFreq[3];
	mQ[ 128] = aK*aOmega*aCodonFreq[6];
	mQ[ 137] = aK*aCodonFreq[15];
	mQ[ 153] = aOmega*aCodonFreq[31];
	mQ[ 169] = aOmega*aCodonFreq[47];
	mQ[ 183] = aOmega*aCodonFreq[0];
	mQ[ 184] = aOmega*aCodonFreq[1];
	mQ[ 185] = aK*aCodonFreq[2];
	mQ[ 186] = -aOmega*aCodonFreq[0]-aOmega*aCodonFreq[1]-aK*aCodonFreq[2]-aK*aOmega*aCodonFreq[7]-aOmega*aCodonFreq[12]-aK*aCodonFreq[16]-aOmega*aCodonFreq[32]-aOmega*aCodonFreq[48];
	mQ[ 190] = aK*aOmega*aCodonFreq[7];
	mQ[ 195] = aOmega*aCodonFreq[12];
	mQ[ 199] = aK*aCodonFreq[16];
	mQ[ 215] = aOmega*aCodonFreq[32];
	mQ[ 231] = aOmega*aCodonFreq[48];
	mQ[ 244] = aK*aOmega*aCodonFreq[0];
	mQ[ 248] = -aK*aOmega*aCodonFreq[0]-aK*aCodonFreq[5]-aCodonFreq[6]-aCodonFreq[7]-aOmega*aCodonFreq[8]-aOmega*aCodonFreq[10]-aK*aOmega*aCodonFreq[17]-aOmega*aCodonFreq[33]-aOmega*aCodonFreq[49];
	mQ[ 249] = aK*aCodonFreq[5];
	mQ[ 250] = aCodonFreq[6];
	mQ[ 251] = aCodonFreq[7];
	mQ[ 252] = aOmega*aCodonFreq[8];
	mQ[ 254] = aOmega*aCodonFreq[10];
	mQ[ 261] = aK*aOmega*aCodonFreq[17];
	mQ[ 277] = aOmega*aCodonFreq[33];
	mQ[ 293] = aOmega*aCodonFreq[49];
	mQ[ 306] = aK*aOmega*aCodonFreq[1];
	mQ[ 309] = aK*aCodonFreq[4];
	mQ[ 310] = -aK*aOmega*aCodonFreq[1]-aK*aCodonFreq[4]-aCodonFreq[6]-aCodonFreq[7]-aOmega*aCodonFreq[9]-aOmega*aCodonFreq[11]-aK*aOmega*aCodonFreq[18]-aOmega*aCodonFreq[34]-aOmega*aCodonFreq[50];
	mQ[ 311] = aCodonFreq[6];
	mQ[ 312] = aCodonFreq[7];
	mQ[ 314] = aOmega*aCodonFreq[9];
	mQ[ 316] = aOmega*aCodonFreq[11];
	mQ[ 323] = aK*aOmega*aCodonFreq[18];
	mQ[ 339] = aOmega*aCodonFreq[34];
	mQ[ 355] = aOmega*aCodonFreq[50];
	mQ[ 368] = aK*aOmega*aCodonFreq[2];
	mQ[ 370] = aCodonFreq[4];
	mQ[ 371] = aCodonFreq[5];
	mQ[ 372] = -aK*aOmega*aCodonFreq[2]-aCodonFreq[4]-aCodonFreq[5]-aK*aCodonFreq[7]-aK*aOmega*aCodonFreq[19]-aOmega*aCodonFreq[35]-aOmega*aCodonFreq[51];
	mQ[ 373] = aK*aCodonFreq[7];
	mQ[ 385] = aK*aOmega*aCodonFreq[19];
	mQ[ 401] = aOmega*aCodonFreq[35];
	mQ[ 417] = aOmega*aCodonFreq[51];
	mQ[ 430] = aK*aOmega*aCodonFreq[3];
	mQ[ 431] = aCodonFreq[4];
	mQ[ 432] = aCodonFreq[5];
	mQ[ 433] = aK*aCodonFreq[6];
	mQ[ 434] = -aK*aOmega*aCodonFreq[3]-aCodonFreq[4]-aCodonFreq[5]-aK*aCodonFreq[6]-aOmega*aCodonFreq[12]-aK*aOmega*aCodonFreq[20]-aOmega*aCodonFreq[36]-aOmega*aCodonFreq[52];
	mQ[ 439] = aOmega*aCodonFreq[12];
	mQ[ 447] = aK*aOmega*aCodonFreq[20];
	mQ[ 463] = aOmega*aCodonFreq[36];
	mQ[ 479] = aOmega*aCodonFreq[52];
	mQ[ 488] = aOmega*aCodonFreq[0];
	mQ[ 492] = aOmega*aCodonFreq[4];
	mQ[ 496] = -aOmega*aCodonFreq[0]-aOmega*aCodonFreq[4]-aK*aCodonFreq[9]-aK*aOmega*aCodonFreq[10]-aK*aOmega*aCodonFreq[21]-aOmega*aCodonFreq[37]-aOmega*aCodonFreq[53];
	mQ[ 497] = aK*aCodonFreq[9];
	mQ[ 498] = aK*aOmega*aCodonFreq[10];
	mQ[ 509] = aK*aOmega*aCodonFreq[21];
	mQ[ 525] = aOmega*aCodonFreq[37];
	mQ[ 541] = aOmega*aCodonFreq[53];
	mQ[ 550] = aOmega*aCodonFreq[1];
	mQ[ 554] = aOmega*aCodonFreq[5];
	mQ[ 557] = aK*aCodonFreq[8];
	mQ[ 558] = -aOmega*aCodonFreq[1]-aOmega*aCodonFreq[5]-aK*aCodonFreq[8]-aK*aOmega*aCodonFreq[11]-aK*aOmega*aCodonFreq[22]-aOmega*aCodonFreq[38]-aOmega*aCodonFreq[54];
	mQ[ 560] = aK*aOmega*aCodonFreq[11];
	mQ[ 571] = aK*aOmega*aCodonFreq[22];
	mQ[ 587] = aOmega*aCodonFreq[38];
	mQ[ 603] = aOmega*aCodonFreq[54];
	mQ[ 610] = aOmega*aCodonFreq[0];
	mQ[ 614] = aOmega*aCodonFreq[4];
	mQ[ 618] = aK*aOmega*aCodonFreq[8];
	mQ[ 620] = -aOmega*aCodonFreq[0]-aOmega*aCodonFreq[4]-aK*aOmega*aCodonFreq[8]-aK*aCodonFreq[11]-aOmega*aCodonFreq[12]-aK*aOmega*aCodonFreq[25]-aOmega*aCodonFreq[41]-aOmega*aCodonFreq[57];
	mQ[ 621] = aK*aCodonFreq[11];
	mQ[ 622] = aOmega*aCodonFreq[12];
	mQ[ 635] = aK*aOmega*aCodonFreq[25];
	mQ[ 651] = aOmega*aCodonFreq[41];
	mQ[ 667] = aOmega*aCodonFreq[57];
	mQ[ 672] = aOmega*aCodonFreq[1];
	mQ[ 676] = aOmega*aCodonFreq[5];
	mQ[ 680] = aK*aOmega*aCodonFreq[9];
	mQ[ 681] = aK*aCodonFreq[10];
	mQ[ 682] = -aOmega*aCodonFreq[1]-aOmega*aCodonFreq[5]-aK*aOmega*aCodonFreq[9]-aK*aCodonFreq[10]-aOmega*aCodonFreq[12]-aK*aOmega*aCodonFreq[26]-aOmega*aCodonFreq[42]-aOmega*aCodonFreq[58];
	mQ[ 683] = aOmega*aCodonFreq[12];
	mQ[ 697] = aK*aOmega*aCodonFreq[26];
	mQ[ 713] = aOmega*aCodonFreq[42];
	mQ[ 729] = aOmega*aCodonFreq[58];
	mQ[ 735] = aOmega*aCodonFreq[3];
	mQ[ 739] = aOmega*aCodonFreq[7];
	mQ[ 742] = aOmega*aCodonFreq[10];
	mQ[ 743] = aOmega*aCodonFreq[11];
	mQ[ 744] = -aOmega*aCodonFreq[3]-aOmega*aCodonFreq[7]-aOmega*aCodonFreq[10]-aOmega*aCodonFreq[11]-aK*aOmega*aCodonFreq[28]-aOmega*aCodonFreq[44]-aOmega*aCodonFreq[60];
	mQ[ 760] = aK*aOmega*aCodonFreq[28];
	mQ[ 776] = aOmega*aCodonFreq[44];
	mQ[ 792] = aOmega*aCodonFreq[60];
	mQ[ 793] = aK*aOmega*aCodonFreq[0];
	mQ[ 806] = -aK*aOmega*aCodonFreq[0]-aK*aCodonFreq[14]-aCodonFreq[15]-aCodonFreq[16]-aK*aOmega*aCodonFreq[17]-aOmega*aCodonFreq[21]-aOmega*aCodonFreq[25]-aOmega*aCodonFreq[29]-aOmega*aCodonFreq[45];
	mQ[ 807] = aK*aCodonFreq[14];
	mQ[ 808] = aCodonFreq[15];
	mQ[ 809] = aCodonFreq[16];
	mQ[ 810] = aK*aOmega*aCodonFreq[17];
	mQ[ 814] = aOmega*aCodonFreq[21];
	mQ[ 818] = aOmega*aCodonFreq[25];
	mQ[ 822] = aOmega*aCodonFreq[29];
	mQ[ 838] = aOmega*aCodonFreq[45];
	mQ[ 855] = aK*aOmega*aCodonFreq[1];
	mQ[ 867] = aK*aCodonFreq[13];
	mQ[ 868] = -aK*aOmega*aCodonFreq[1]-aK*aCodonFreq[13]-aCodonFreq[15]-aCodonFreq[16]-aK*aOmega*aCodonFreq[18]-aOmega*aCodonFreq[22]-aOmega*aCodonFreq[26]-aOmega*aCodonFreq[30]-aOmega*aCodonFreq[46];
	mQ[ 869] = aCodonFreq[15];
	mQ[ 870] = aCodonFreq[16];
	mQ[ 872] = aK*aOmega*aCodonFreq[18];
	mQ[ 876] = aOmega*aCodonFreq[22];
	mQ[ 880] = aOmega*aCodonFreq[26];
	mQ[ 884] = aOmega*aCodonFreq[30];
	mQ[ 900] = aOmega*aCodonFreq[46];
	mQ[ 917] = aK*aCodonFreq[2];
	mQ[ 928] = aCodonFreq[13];
	mQ[ 929] = aCodonFreq[14];
	mQ[ 930] = -aK*aCodonFreq[2]-aCodonFreq[13]-aCodonFreq[14]-aK*aCodonFreq[16]-aK*aOmega*aCodonFreq[19]-aOmega*aCodonFreq[23]-aOmega*aCodonFreq[27]-aOmega*aCodonFreq[31]-aOmega*aCodonFreq[47];
	mQ[ 931] = aK*aCodonFreq[16];
	mQ[ 934] = aK*aOmega*aCodonFreq[19];
	mQ[ 938] = aOmega*aCodonFreq[23];
	mQ[ 942] = aOmega*aCodonFreq[27];
	mQ[ 946] = aOmega*aCodonFreq[31];
	mQ[ 962] = aOmega*aCodonFreq[47];
	mQ[ 979] = aK*aCodonFreq[3];
	mQ[ 989] = aCodonFreq[13];
	mQ[ 990] = aCodonFreq[14];
	mQ[ 991] = aK*aCodonFreq[15];
	mQ[ 992] = -aK*aCodonFreq[3]-aCodonFreq[13]-aCodonFreq[14]-aK*aCodonFreq[15]-aK*aOmega*aCodonFreq[20]-aOmega*aCodonFreq[24]-aOmega*aCodonFreq[28]-aOmega*aCodonFreq[32]-aOmega*aCodonFreq[48];
	mQ[ 996] = aK*aOmega*aCodonFreq[20];
	mQ[1000] = aOmega*aCodonFreq[24];
	mQ[1004] = aOmega*aCodonFreq[28];
	mQ[1008] = aOmega*aCodonFreq[32];
	mQ[1024] = aOmega*aCodonFreq[48];
	mQ[1041] = aK*aOmega*aCodonFreq[4];
	mQ[1050] = aK*aOmega*aCodonFreq[13];
	mQ[1054] = -aK*aOmega*aCodonFreq[4]-aK*aOmega*aCodonFreq[13]-aK*aCodonFreq[18]-aCodonFreq[19]-aCodonFreq[20]-aOmega*aCodonFreq[21]-aOmega*aCodonFreq[25]-aOmega*aCodonFreq[33]-aOmega*aCodonFreq[49];
	mQ[1055] = aK*aCodonFreq[18];
	mQ[1056] = aCodonFreq[19];
	mQ[1057] = aCodonFreq[20];
	mQ[1058] = aOmega*aCodonFreq[21];
	mQ[1062] = aOmega*aCodonFreq[25];
	mQ[1070] = aOmega*aCodonFreq[33];
	mQ[1086] = aOmega*aCodonFreq[49];
	mQ[1103] = aK*aOmega*aCodonFreq[5];
	mQ[1112] = aK*aOmega*aCodonFreq[14];
	mQ[1115] = aK*aCodonFreq[17];
	mQ[1116] = -aK*aOmega*aCodonFreq[5]-aK*aOmega*aCodonFreq[14]-aK*aCodonFreq[17]-aCodonFreq[19]-aCodonFreq[20]-aOmega*aCodonFreq[22]-aOmega*aCodonFreq[26]-aOmega*aCodonFreq[34]-aOmega*aCodonFreq[50];
	mQ[1117] = aCodonFreq[19];
	mQ[1118] = aCodonFreq[20];
	mQ[1120] = aOmega*aCodonFreq[22];
	mQ[1124] = aOmega*aCodonFreq[26];
	mQ[1132] = aOmega*aCodonFreq[34];
	mQ[1148] = aOmega*aCodonFreq[50];
	mQ[1165] = aK*aOmega*aCodonFreq[6];
	mQ[1174] = aK*aOmega*aCodonFreq[15];
	mQ[1176] = aCodonFreq[17];
	mQ[1177] = aCodonFreq[18];
	mQ[1178] = -aK*aOmega*aCodonFreq[6]-aK*aOmega*aCodonFreq[15]-aCodonFreq[17]-aCodonFreq[18]-aK*aCodonFreq[20]-aOmega*aCodonFreq[23]-aOmega*aCodonFreq[27]-aOmega*aCodonFreq[35]-aOmega*aCodonFreq[51];
	mQ[1179] = aK*aCodonFreq[20];
	mQ[1182] = aOmega*aCodonFreq[23];
	mQ[1186] = aOmega*aCodonFreq[27];
	mQ[1194] = aOmega*aCodonFreq[35];
	mQ[1210] = aOmega*aCodonFreq[51];
	mQ[1227] = aK*aOmega*aCodonFreq[7];
	mQ[1236] = aK*aOmega*aCodonFreq[16];
	mQ[1237] = aCodonFreq[17];
	mQ[1238] = aCodonFreq[18];
	mQ[1239] = aK*aCodonFreq[19];
	mQ[1240] = -aK*aOmega*aCodonFreq[7]-aK*aOmega*aCodonFreq[16]-aCodonFreq[17]-aCodonFreq[18]-aK*aCodonFreq[19]-aOmega*aCodonFreq[24]-aOmega*aCodonFreq[28]-aOmega*aCodonFreq[36]-aOmega*aCodonFreq[52];
	mQ[1244] = aOmega*aCodonFreq[24];
	mQ[1248] = aOmega*aCodonFreq[28];
	mQ[1256] = aOmega*aCodonFreq[36];
	mQ[1272] = aOmega*aCodonFreq[52];
	mQ[1289] = aK*aOmega*aCodonFreq[8];
	mQ[1294] = aOmega*aCodonFreq[13];
	mQ[1298] = aOmega*aCodonFreq[17];
	mQ[1302] = -aK*aOmega*aCodonFreq[8]-aOmega*aCodonFreq[13]-aOmega*aCodonFreq[17]-aK*aCodonFreq[22]-aOmega*aCodonFreq[23]-aOmega*aCodonFreq[24]-aK*aOmega*aCodonFreq[25]-aOmega*aCodonFreq[37]-aOmega*aCodonFreq[53];
	mQ[1303] = aK*aCodonFreq[22];
	mQ[1304] = aOmega*aCodonFreq[23];
	mQ[1305] = aOmega*aCodonFreq[24];
	mQ[1306] = aK*aOmega*aCodonFreq[25];
	mQ[1318] = aOmega*aCodonFreq[37];
	mQ[1334] = aOmega*aCodonFreq[53];
	mQ[1351] = aK*aOmega*aCodonFreq[9];
	mQ[1356] = aOmega*aCodonFreq[14];
	mQ[1360] = aOmega*aCodonFreq[18];
	mQ[1363] = aK*aCodonFreq[21];
	mQ[1364] = -aK*aOmega*aCodonFreq[9]-aOmega*aCodonFreq[14]-aOmega*aCodonFreq[18]-aK*aCodonFreq[21]-aOmega*aCodonFreq[23]-aOmega*aCodonFreq[24]-aK*aOmega*aCodonFreq[26]-aOmega*aCodonFreq[38]-aOmega*aCodonFreq[54];
	mQ[1365] = aOmega*aCodonFreq[23];
	mQ[1366] = aOmega*aCodonFreq[24];
	mQ[1368] = aK*aOmega*aCodonFreq[26];
	mQ[1380] = aOmega*aCodonFreq[38];
	mQ[1396] = aOmega*aCodonFreq[54];
	mQ[1418] = aOmega*aCodonFreq[15];
	mQ[1422] = aOmega*aCodonFreq[19];
	mQ[1424] = aOmega*aCodonFreq[21];
	mQ[1425] = aOmega*aCodonFreq[22];
	mQ[1426] = -aOmega*aCodonFreq[15]-aOmega*aCodonFreq[19]-aOmega*aCodonFreq[21]-aOmega*aCodonFreq[22]-aK*aCodonFreq[24]-aK*aOmega*aCodonFreq[27]-aOmega*aCodonFreq[39]-aOmega*aCodonFreq[55];
	mQ[1427] = aK*aCodonFreq[24];
	mQ[1430] = aK*aOmega*aCodonFreq[27];
	mQ[1442] = aOmega*aCodonFreq[39];
	mQ[1458] = aOmega*aCodonFreq[55];
	mQ[1480] = aOmega*aCodonFreq[16];
	mQ[1484] = aOmega*aCodonFreq[20];
	mQ[1485] = aOmega*aCodonFreq[21];
	mQ[1486] = aOmega*aCodonFreq[22];
	mQ[1487] = aK*aCodonFreq[23];
	mQ[1488] = -aOmega*aCodonFreq[16]-aOmega*aCodonFreq[20]-aOmega*aCodonFreq[21]-aOmega*aCodonFreq[22]-aK*aCodonFreq[23]-aK*aOmega*aCodonFreq[28]-aOmega*aCodonFreq[40]-aOmega*aCodonFreq[56];
	mQ[1492] = aK*aOmega*aCodonFreq[28];
	mQ[1504] = aOmega*aCodonFreq[40];
	mQ[1520] = aOmega*aCodonFreq[56];
	mQ[1535] = aK*aOmega*aCodonFreq[10];
	mQ[1538] = aOmega*aCodonFreq[13];
	mQ[1542] = aOmega*aCodonFreq[17];
	mQ[1546] = aK*aOmega*aCodonFreq[21];
	mQ[1550] = -aK*aOmega*aCodonFreq[10]-aOmega*aCodonFreq[13]-aOmega*aCodonFreq[17]-aK*aOmega*aCodonFreq[21]-aK*aCodonFreq[26]-aCodonFreq[27]-aCodonFreq[28]-aOmega*aCodonFreq[41]-aOmega*aCodonFreq[57];
	mQ[1551] = aK*aCodonFreq[26];
	mQ[1552] = aCodonFreq[27];
	mQ[1553] = aCodonFreq[28];
	mQ[1566] = aOmega*aCodonFreq[41];
	mQ[1582] = aOmega*aCodonFreq[57];
	mQ[1597] = aK*aOmega*aCodonFreq[11];
	mQ[1600] = aOmega*aCodonFreq[14];
	mQ[1604] = aOmega*aCodonFreq[18];
	mQ[1608] = aK*aOmega*aCodonFreq[22];
	mQ[1611] = aK*aCodonFreq[25];
	mQ[1612] = -aK*aOmega*aCodonFreq[11]-aOmega*aCodonFreq[14]-aOmega*aCodonFreq[18]-aK*aOmega*aCodonFreq[22]-aK*aCodonFreq[25]-aCodonFreq[27]-aCodonFreq[28]-aOmega*aCodonFreq[42]-aOmega*aCodonFreq[58];
	mQ[1613] = aCodonFreq[27];
	mQ[1614] = aCodonFreq[28];
	mQ[1628] = aOmega*aCodonFreq[42];
	mQ[1644] = aOmega*aCodonFreq[58];
	mQ[1662] = aOmega*aCodonFreq[15];
	mQ[1666] = aOmega*aCodonFreq[19];
	mQ[1670] = aK*aOmega*aCodonFreq[23];
	mQ[1672] = aCodonFreq[25];
	mQ[1673] = aCodonFreq[26];
	mQ[1674] = -aOmega*aCodonFreq[15]-aOmega*aCodonFreq[19]-aK*aOmega*aCodonFreq[23]-aCodonFreq[25]-aCodonFreq[26]-aK*aCodonFreq[28]-aCodonFreq[43]-aOmega*aCodonFreq[59];
	mQ[1675] = aK*aCodonFreq[28];
	mQ[1690] = aCodonFreq[43];
	mQ[1706] = aOmega*aCodonFreq[59];
	mQ[1720] = aK*aOmega*aCodonFreq[12];
	mQ[1724] = aOmega*aCodonFreq[16];
	mQ[1728] = aOmega*aCodonFreq[20];
	mQ[1732] = aK*aOmega*aCodonFreq[24];
	mQ[1733] = aCodonFreq[25];
	mQ[1734] = aCodonFreq[26];
	mQ[1735] = aK*aCodonFreq[27];
	mQ[1736] = -aK*aOmega*aCodonFreq[12]-aOmega*aCodonFreq[16]-aOmega*aCodonFreq[20]-aK*aOmega*aCodonFreq[24]-aCodonFreq[25]-aCodonFreq[26]-aK*aCodonFreq[27]-aCodonFreq[44]-aOmega*aCodonFreq[60];
	mQ[1752] = aCodonFreq[44];
	mQ[1768] = aOmega*aCodonFreq[60];
	mQ[1769] = aOmega*aCodonFreq[0];
	mQ[1782] = aOmega*aCodonFreq[13];
	mQ[1798] = -aOmega*aCodonFreq[0]-aOmega*aCodonFreq[13]-aK*aCodonFreq[30]-aCodonFreq[31]-aOmega*aCodonFreq[32]-aK*aOmega*aCodonFreq[33]-aOmega*aCodonFreq[37]-aOmega*aCodonFreq[41]-aK*aOmega*aCodonFreq[45];
	mQ[1799] = aK*aCodonFreq[30];
	mQ[1800] = aCodonFreq[31];
	mQ[1801] = aOmega*aCodonFreq[32];
	mQ[1802] = aK*aOmega*aCodonFreq[33];
	mQ[1806] = aOmega*aCodonFreq[37];
	mQ[1810] = aOmega*aCodonFreq[41];
	mQ[1814] = aK*aOmega*aCodonFreq[45];
	mQ[1831] = aOmega*aCodonFreq[1];
	mQ[1844] = aOmega*aCodonFreq[14];
	mQ[1859] = aK*aCodonFreq[29];
	mQ[1860] = -aOmega*aCodonFreq[1]-aOmega*aCodonFreq[14]-aK*aCodonFreq[29]-aCodonFreq[31]-aOmega*aCodonFreq[32]-aK*aOmega*aCodonFreq[34]-aOmega*aCodonFreq[38]-aOmega*aCodonFreq[42]-aK*aOmega*aCodonFreq[46];
	mQ[1861] = aCodonFreq[31];
	mQ[1862] = aOmega*aCodonFreq[32];
	mQ[1864] = aK*aOmega*aCodonFreq[34];
	mQ[1868] = aOmega*aCodonFreq[38];
	mQ[1872] = aOmega*aCodonFreq[42];
	mQ[1876] = aK*aOmega*aCodonFreq[46];
	mQ[1893] = aOmega*aCodonFreq[2];
	mQ[1906] = aOmega*aCodonFreq[15];
	mQ[1920] = aCodonFreq[29];
	mQ[1921] = aCodonFreq[30];
	mQ[1922] = -aOmega*aCodonFreq[2]-aOmega*aCodonFreq[15]-aCodonFreq[29]-aCodonFreq[30]-aK*aOmega*aCodonFreq[32]-aK*aOmega*aCodonFreq[35]-aOmega*aCodonFreq[39]-aOmega*aCodonFreq[43]-aK*aOmega*aCodonFreq[47];
	mQ[1923] = aK*aOmega*aCodonFreq[32];
	mQ[1926] = aK*aOmega*aCodonFreq[35];
	mQ[1930] = aOmega*aCodonFreq[39];
	mQ[1934] = aOmega*aCodonFreq[43];
	mQ[1938] = aK*aOmega*aCodonFreq[47];
	mQ[1955] = aOmega*aCodonFreq[3];
	mQ[1968] = aOmega*aCodonFreq[16];
	mQ[1981] = aOmega*aCodonFreq[29];
	mQ[1982] = aOmega*aCodonFreq[30];
	mQ[1983] = aK*aOmega*aCodonFreq[31];
	mQ[1984] = -aOmega*aCodonFreq[3]-aOmega*aCodonFreq[16]-aOmega*aCodonFreq[29]-aOmega*aCodonFreq[30]-aK*aOmega*aCodonFreq[31]-aK*aOmega*aCodonFreq[36]-aOmega*aCodonFreq[40]-aOmega*aCodonFreq[44]-aK*aOmega*aCodonFreq[48];
	mQ[1988] = aK*aOmega*aCodonFreq[36];
	mQ[1992] = aOmega*aCodonFreq[40];
	mQ[1996] = aOmega*aCodonFreq[44];
	mQ[2000] = aK*aOmega*aCodonFreq[48];
	mQ[2017] = aOmega*aCodonFreq[4];
	mQ[2030] = aOmega*aCodonFreq[17];
	mQ[2042] = aK*aOmega*aCodonFreq[29];
	mQ[2046] = -aOmega*aCodonFreq[4]-aOmega*aCodonFreq[17]-aK*aOmega*aCodonFreq[29]-aK*aCodonFreq[34]-aCodonFreq[35]-aCodonFreq[36]-aOmega*aCodonFreq[37]-aOmega*aCodonFreq[41]-aK*aOmega*aCodonFreq[49];
	mQ[2047] = aK*aCodonFreq[34];
	mQ[2048] = aCodonFreq[35];
	mQ[2049] = aCodonFreq[36];
	mQ[2050] = aOmega*aCodonFreq[37];
	mQ[2054] = aOmega*aCodonFreq[41];
	mQ[2062] = aK*aOmega*aCodonFreq[49];
	mQ[2079] = aOmega*aCodonFreq[5];
	mQ[2092] = aOmega*aCodonFreq[18];
	mQ[2104] = aK*aOmega*aCodonFreq[30];
	mQ[2107] = aK*aCodonFreq[33];
	mQ[2108] = -aOmega*aCodonFreq[5]-aOmega*aCodonFreq[18]-aK*aOmega*aCodonFreq[30]-aK*aCodonFreq[33]-aCodonFreq[35]-aCodonFreq[36]-aOmega*aCodonFreq[38]-aOmega*aCodonFreq[42]-aK*aOmega*aCodonFreq[50];
	mQ[2109] = aCodonFreq[35];
	mQ[2110] = aCodonFreq[36];
	mQ[2112] = aOmega*aCodonFreq[38];
	mQ[2116] = aOmega*aCodonFreq[42];
	mQ[2124] = aK*aOmega*aCodonFreq[50];
	mQ[2141] = aOmega*aCodonFreq[6];
	mQ[2154] = aOmega*aCodonFreq[19];
	mQ[2166] = aK*aOmega*aCodonFreq[31];
	mQ[2168] = aCodonFreq[33];
	mQ[2169] = aCodonFreq[34];
	mQ[2170] = -aOmega*aCodonFreq[6]-aOmega*aCodonFreq[19]-aK*aOmega*aCodonFreq[31]-aCodonFreq[33]-aCodonFreq[34]-aK*aCodonFreq[36]-aOmega*aCodonFreq[39]-aOmega*aCodonFreq[43]-aK*aOmega*aCodonFreq[51];
	mQ[2171] = aK*aCodonFreq[36];
	mQ[2174] = aOmega*aCodonFreq[39];
	mQ[2178] = aOmega*aCodonFreq[43];
	mQ[2186] = aK*aOmega*aCodonFreq[51];
	mQ[2203] = aOmega*aCodonFreq[7];
	mQ[2216] = aOmega*aCodonFreq[20];
	mQ[2228] = aK*aOmega*aCodonFreq[32];
	mQ[2229] = aCodonFreq[33];
	mQ[2230] = aCodonFreq[34];
	mQ[2231] = aK*aCodonFreq[35];
	mQ[2232] = -aOmega*aCodonFreq[7]-aOmega*aCodonFreq[20]-aK*aOmega*aCodonFreq[32]-aCodonFreq[33]-aCodonFreq[34]-aK*aCodonFreq[35]-aOmega*aCodonFreq[40]-aOmega*aCodonFreq[44]-aK*aOmega*aCodonFreq[52];
	mQ[2236] = aOmega*aCodonFreq[40];
	mQ[2240] = aOmega*aCodonFreq[44];
	mQ[2248] = aK*aOmega*aCodonFreq[52];
	mQ[2265] = aOmega*aCodonFreq[8];
	mQ[2278] = aOmega*aCodonFreq[21];
	mQ[2286] = aOmega*aCodonFreq[29];
	mQ[2290] = aOmega*aCodonFreq[33];
	mQ[2294] = -aOmega*aCodonFreq[8]-aOmega*aCodonFreq[21]-aOmega*aCodonFreq[29]-aOmega*aCodonFreq[33]-aK*aCodonFreq[38]-aOmega*aCodonFreq[39]-aOmega*aCodonFreq[40]-aK*aOmega*aCodonFreq[41]-aK*aOmega*aCodonFreq[53];
	mQ[2295] = aK*aCodonFreq[38];
	mQ[2296] = aOmega*aCodonFreq[39];
	mQ[2297] = aOmega*aCodonFreq[40];
	mQ[2298] = aK*aOmega*aCodonFreq[41];
	mQ[2310] = aK*aOmega*aCodonFreq[53];
	mQ[2327] = aOmega*aCodonFreq[9];
	mQ[2340] = aOmega*aCodonFreq[22];
	mQ[2348] = aOmega*aCodonFreq[30];
	mQ[2352] = aOmega*aCodonFreq[34];
	mQ[2355] = aK*aCodonFreq[37];
	mQ[2356] = -aOmega*aCodonFreq[9]-aOmega*aCodonFreq[22]-aOmega*aCodonFreq[30]-aOmega*aCodonFreq[34]-aK*aCodonFreq[37]-aOmega*aCodonFreq[39]-aOmega*aCodonFreq[40]-aK*aOmega*aCodonFreq[42]-aK*aOmega*aCodonFreq[54];
	mQ[2357] = aOmega*aCodonFreq[39];
	mQ[2358] = aOmega*aCodonFreq[40];
	mQ[2360] = aK*aOmega*aCodonFreq[42];
	mQ[2372] = aK*aOmega*aCodonFreq[54];
	mQ[2402] = aOmega*aCodonFreq[23];
	mQ[2410] = aOmega*aCodonFreq[31];
	mQ[2414] = aOmega*aCodonFreq[35];
	mQ[2416] = aOmega*aCodonFreq[37];
	mQ[2417] = aOmega*aCodonFreq[38];
	mQ[2418] = -aOmega*aCodonFreq[23]-aOmega*aCodonFreq[31]-aOmega*aCodonFreq[35]-aOmega*aCodonFreq[37]-aOmega*aCodonFreq[38]-aK*aCodonFreq[40]-aK*aOmega*aCodonFreq[43]-aK*aOmega*aCodonFreq[55];
	mQ[2419] = aK*aCodonFreq[40];
	mQ[2422] = aK*aOmega*aCodonFreq[43];
	mQ[2434] = aK*aOmega*aCodonFreq[55];
	mQ[2464] = aOmega*aCodonFreq[24];
	mQ[2472] = aOmega*aCodonFreq[32];
	mQ[2476] = aOmega*aCodonFreq[36];
	mQ[2477] = aOmega*aCodonFreq[37];
	mQ[2478] = aOmega*aCodonFreq[38];
	mQ[2479] = aK*aCodonFreq[39];
	mQ[2480] = -aOmega*aCodonFreq[24]-aOmega*aCodonFreq[32]-aOmega*aCodonFreq[36]-aOmega*aCodonFreq[37]-aOmega*aCodonFreq[38]-aK*aCodonFreq[39]-aK*aOmega*aCodonFreq[44]-aK*aOmega*aCodonFreq[56];
	mQ[2484] = aK*aOmega*aCodonFreq[44];
	mQ[2496] = aK*aOmega*aCodonFreq[56];
	mQ[2511] = aOmega*aCodonFreq[10];
	mQ[2526] = aOmega*aCodonFreq[25];
	mQ[2530] = aOmega*aCodonFreq[29];
	mQ[2534] = aOmega*aCodonFreq[33];
	mQ[2538] = aK*aOmega*aCodonFreq[37];
	mQ[2542] = -aOmega*aCodonFreq[10]-aOmega*aCodonFreq[25]-aOmega*aCodonFreq[29]-aOmega*aCodonFreq[33]-aK*aOmega*aCodonFreq[37]-aK*aCodonFreq[42]-aOmega*aCodonFreq[43]-aOmega*aCodonFreq[44]-aK*aOmega*aCodonFreq[57];
	mQ[2543] = aK*aCodonFreq[42];
	mQ[2544] = aOmega*aCodonFreq[43];
	mQ[2545] = aOmega*aCodonFreq[44];
	mQ[2558] = aK*aOmega*aCodonFreq[57];
	mQ[2573] = aOmega*aCodonFreq[11];
	mQ[2588] = aOmega*aCodonFreq[26];
	mQ[2592] = aOmega*aCodonFreq[30];
	mQ[2596] = aOmega*aCodonFreq[34];
	mQ[2600] = aK*aOmega*aCodonFreq[38];
	mQ[2603] = aK*aCodonFreq[41];
	mQ[2604] = -aOmega*aCodonFreq[11]-aOmega*aCodonFreq[26]-aOmega*aCodonFreq[30]-aOmega*aCodonFreq[34]-aK*aOmega*aCodonFreq[38]-aK*aCodonFreq[41]-aOmega*aCodonFreq[43]-aOmega*aCodonFreq[44]-aK*aOmega*aCodonFreq[58];
	mQ[2605] = aOmega*aCodonFreq[43];
	mQ[2606] = aOmega*aCodonFreq[44];
	mQ[2620] = aK*aOmega*aCodonFreq[58];
	mQ[2650] = aCodonFreq[27];
	mQ[2654] = aOmega*aCodonFreq[31];
	mQ[2658] = aOmega*aCodonFreq[35];
	mQ[2662] = aK*aOmega*aCodonFreq[39];
	mQ[2664] = aOmega*aCodonFreq[41];
	mQ[2665] = aOmega*aCodonFreq[42];
	mQ[2666] = -aCodonFreq[27]-aOmega*aCodonFreq[31]-aOmega*aCodonFreq[35]-aK*aOmega*aCodonFreq[39]-aOmega*aCodonFreq[41]-aOmega*aCodonFreq[42]-aK*aCodonFreq[44]-aK*aOmega*aCodonFreq[59];
	mQ[2667] = aK*aCodonFreq[44];
	mQ[2682] = aK*aOmega*aCodonFreq[59];
	mQ[2696] = aOmega*aCodonFreq[12];
	mQ[2712] = aCodonFreq[28];
	mQ[2716] = aOmega*aCodonFreq[32];
	mQ[2720] = aOmega*aCodonFreq[36];
	mQ[2724] = aK*aOmega*aCodonFreq[40];
	mQ[2725] = aOmega*aCodonFreq[41];
	mQ[2726] = aOmega*aCodonFreq[42];
	mQ[2727] = aK*aCodonFreq[43];
	mQ[2728] = -aOmega*aCodonFreq[12]-aCodonFreq[28]-aOmega*aCodonFreq[32]-aOmega*aCodonFreq[36]-aK*aOmega*aCodonFreq[40]-aOmega*aCodonFreq[41]-aOmega*aCodonFreq[42]-aK*aCodonFreq[43]-aK*aOmega*aCodonFreq[60];
	mQ[2744] = aK*aOmega*aCodonFreq[60];
	mQ[2745] = aOmega*aCodonFreq[0];
	mQ[2758] = aOmega*aCodonFreq[13];
	mQ[2774] = aK*aOmega*aCodonFreq[29];
	mQ[2790] = -aOmega*aCodonFreq[0]-aOmega*aCodonFreq[13]-aK*aOmega*aCodonFreq[29]-aK*aCodonFreq[46]-aCodonFreq[47]-aCodonFreq[48]-aK*aOmega*aCodonFreq[49]-aOmega*aCodonFreq[53]-aOmega*aCodonFreq[57];
	mQ[2791] = aK*aCodonFreq[46];
	mQ[2792] = aCodonFreq[47];
	mQ[2793] = aCodonFreq[48];
	mQ[2794] = aK*aOmega*aCodonFreq[49];
	mQ[2798] = aOmega*aCodonFreq[53];
	mQ[2802] = aOmega*aCodonFreq[57];
	mQ[2807] = aOmega*aCodonFreq[1];
	mQ[2820] = aOmega*aCodonFreq[14];
	mQ[2836] = aK*aOmega*aCodonFreq[30];
	mQ[2851] = aK*aCodonFreq[45];
	mQ[2852] = -aOmega*aCodonFreq[1]-aOmega*aCodonFreq[14]-aK*aOmega*aCodonFreq[30]-aK*aCodonFreq[45]-aCodonFreq[47]-aCodonFreq[48]-aK*aOmega*aCodonFreq[50]-aOmega*aCodonFreq[54]-aOmega*aCodonFreq[58];
	mQ[2853] = aCodonFreq[47];
	mQ[2854] = aCodonFreq[48];
	mQ[2856] = aK*aOmega*aCodonFreq[50];
	mQ[2860] = aOmega*aCodonFreq[54];
	mQ[2864] = aOmega*aCodonFreq[58];
	mQ[2869] = aOmega*aCodonFreq[2];
	mQ[2882] = aOmega*aCodonFreq[15];
	mQ[2898] = aK*aOmega*aCodonFreq[31];
	mQ[2912] = aCodonFreq[45];
	mQ[2913] = aCodonFreq[46];
	mQ[2914] = -aOmega*aCodonFreq[2]-aOmega*aCodonFreq[15]-aK*aOmega*aCodonFreq[31]-aCodonFreq[45]-aCodonFreq[46]-aK*aCodonFreq[48]-aK*aOmega*aCodonFreq[51]-aOmega*aCodonFreq[55]-aOmega*aCodonFreq[59];
	mQ[2915] = aK*aCodonFreq[48];
	mQ[2918] = aK*aOmega*aCodonFreq[51];
	mQ[2922] = aOmega*aCodonFreq[55];
	mQ[2926] = aOmega*aCodonFreq[59];
	mQ[2931] = aOmega*aCodonFreq[3];
	mQ[2944] = aOmega*aCodonFreq[16];
	mQ[2960] = aK*aOmega*aCodonFreq[32];
	mQ[2973] = aCodonFreq[45];
	mQ[2974] = aCodonFreq[46];
	mQ[2975] = aK*aCodonFreq[47];
	mQ[2976] = -aOmega*aCodonFreq[3]-aOmega*aCodonFreq[16]-aK*aOmega*aCodonFreq[32]-aCodonFreq[45]-aCodonFreq[46]-aK*aCodonFreq[47]-aK*aOmega*aCodonFreq[52]-aOmega*aCodonFreq[56]-aOmega*aCodonFreq[60];
	mQ[2980] = aK*aOmega*aCodonFreq[52];
	mQ[2984] = aOmega*aCodonFreq[56];
	mQ[2988] = aOmega*aCodonFreq[60];
	mQ[2993] = aOmega*aCodonFreq[4];
	mQ[3006] = aOmega*aCodonFreq[17];
	mQ[3022] = aK*aOmega*aCodonFreq[33];
	mQ[3034] = aK*aOmega*aCodonFreq[45];
	mQ[3038] = -aOmega*aCodonFreq[4]-aOmega*aCodonFreq[17]-aK*aOmega*aCodonFreq[33]-aK*aOmega*aCodonFreq[45]-aK*aCodonFreq[50]-aCodonFreq[51]-aCodonFreq[52]-aOmega*aCodonFreq[53]-aOmega*aCodonFreq[57];
	mQ[3039] = aK*aCodonFreq[50];
	mQ[3040] = aCodonFreq[51];
	mQ[3041] = aCodonFreq[52];
	mQ[3042] = aOmega*aCodonFreq[53];
	mQ[3046] = aOmega*aCodonFreq[57];
	mQ[3055] = aOmega*aCodonFreq[5];
	mQ[3068] = aOmega*aCodonFreq[18];
	mQ[3084] = aK*aOmega*aCodonFreq[34];
	mQ[3096] = aK*aOmega*aCodonFreq[46];
	mQ[3099] = aK*aCodonFreq[49];
	mQ[3100] = -aOmega*aCodonFreq[5]-aOmega*aCodonFreq[18]-aK*aOmega*aCodonFreq[34]-aK*aOmega*aCodonFreq[46]-aK*aCodonFreq[49]-aCodonFreq[51]-aCodonFreq[52]-aOmega*aCodonFreq[54]-aOmega*aCodonFreq[58];
	mQ[3101] = aCodonFreq[51];
	mQ[3102] = aCodonFreq[52];
	mQ[3104] = aOmega*aCodonFreq[54];
	mQ[3108] = aOmega*aCodonFreq[58];
	mQ[3117] = aOmega*aCodonFreq[6];
	mQ[3130] = aOmega*aCodonFreq[19];
	mQ[3146] = aK*aOmega*aCodonFreq[35];
	mQ[3158] = aK*aOmega*aCodonFreq[47];
	mQ[3160] = aCodonFreq[49];
	mQ[3161] = aCodonFreq[50];
	mQ[3162] = -aOmega*aCodonFreq[6]-aOmega*aCodonFreq[19]-aK*aOmega*aCodonFreq[35]-aK*aOmega*aCodonFreq[47]-aCodonFreq[49]-aCodonFreq[50]-aK*aCodonFreq[52]-aOmega*aCodonFreq[55]-aOmega*aCodonFreq[59];
	mQ[3163] = aK*aCodonFreq[52];
	mQ[3166] = aOmega*aCodonFreq[55];
	mQ[3170] = aOmega*aCodonFreq[59];
	mQ[3179] = aOmega*aCodonFreq[7];
	mQ[3192] = aOmega*aCodonFreq[20];
	mQ[3208] = aK*aOmega*aCodonFreq[36];
	mQ[3220] = aK*aOmega*aCodonFreq[48];
	mQ[3221] = aCodonFreq[49];
	mQ[3222] = aCodonFreq[50];
	mQ[3223] = aK*aCodonFreq[51];
	mQ[3224] = -aOmega*aCodonFreq[7]-aOmega*aCodonFreq[20]-aK*aOmega*aCodonFreq[36]-aK*aOmega*aCodonFreq[48]-aCodonFreq[49]-aCodonFreq[50]-aK*aCodonFreq[51]-aOmega*aCodonFreq[56]-aOmega*aCodonFreq[60];
	mQ[3228] = aOmega*aCodonFreq[56];
	mQ[3232] = aOmega*aCodonFreq[60];
	mQ[3241] = aOmega*aCodonFreq[8];
	mQ[3254] = aOmega*aCodonFreq[21];
	mQ[3270] = aK*aOmega*aCodonFreq[37];
	mQ[3278] = aOmega*aCodonFreq[45];
	mQ[3282] = aOmega*aCodonFreq[49];
	mQ[3286] = -aOmega*aCodonFreq[8]-aOmega*aCodonFreq[21]-aK*aOmega*aCodonFreq[37]-aOmega*aCodonFreq[45]-aOmega*aCodonFreq[49]-aK*aCodonFreq[54]-aOmega*aCodonFreq[55]-aOmega*aCodonFreq[56]-aK*aOmega*aCodonFreq[57];
	mQ[3287] = aK*aCodonFreq[54];
	mQ[3288] = aOmega*aCodonFreq[55];
	mQ[3289] = aOmega*aCodonFreq[56];
	mQ[3290] = aK*aOmega*aCodonFreq[57];
	mQ[3303] = aOmega*aCodonFreq[9];
	mQ[3316] = aOmega*aCodonFreq[22];
	mQ[3332] = aK*aOmega*aCodonFreq[38];
	mQ[3340] = aOmega*aCodonFreq[46];
	mQ[3344] = aOmega*aCodonFreq[50];
	mQ[3347] = aK*aCodonFreq[53];
	mQ[3348] = -aOmega*aCodonFreq[9]-aOmega*aCodonFreq[22]-aK*aOmega*aCodonFreq[38]-aOmega*aCodonFreq[46]-aOmega*aCodonFreq[50]-aK*aCodonFreq[53]-aOmega*aCodonFreq[55]-aOmega*aCodonFreq[56]-aK*aOmega*aCodonFreq[58];
	mQ[3349] = aOmega*aCodonFreq[55];
	mQ[3350] = aOmega*aCodonFreq[56];
	mQ[3352] = aK*aOmega*aCodonFreq[58];
	mQ[3378] = aOmega*aCodonFreq[23];
	mQ[3394] = aK*aOmega*aCodonFreq[39];
	mQ[3402] = aOmega*aCodonFreq[47];
	mQ[3406] = aOmega*aCodonFreq[51];
	mQ[3408] = aOmega*aCodonFreq[53];
	mQ[3409] = aOmega*aCodonFreq[54];
	mQ[3410] = -aOmega*aCodonFreq[23]-aK*aOmega*aCodonFreq[39]-aOmega*aCodonFreq[47]-aOmega*aCodonFreq[51]-aOmega*aCodonFreq[53]-aOmega*aCodonFreq[54]-aK*aCodonFreq[56]-aK*aOmega*aCodonFreq[59];
	mQ[3411] = aK*aCodonFreq[56];
	mQ[3414] = aK*aOmega*aCodonFreq[59];
	mQ[3440] = aOmega*aCodonFreq[24];
	mQ[3456] = aK*aOmega*aCodonFreq[40];
	mQ[3464] = aOmega*aCodonFreq[48];
	mQ[3468] = aOmega*aCodonFreq[52];
	mQ[3469] = aOmega*aCodonFreq[53];
	mQ[3470] = aOmega*aCodonFreq[54];
	mQ[3471] = aK*aCodonFreq[55];
	mQ[3472] = -aOmega*aCodonFreq[24]-aK*aOmega*aCodonFreq[40]-aOmega*aCodonFreq[48]-aOmega*aCodonFreq[52]-aOmega*aCodonFreq[53]-aOmega*aCodonFreq[54]-aK*aCodonFreq[55]-aK*aOmega*aCodonFreq[60];
	mQ[3476] = aK*aOmega*aCodonFreq[60];
	mQ[3487] = aOmega*aCodonFreq[10];
	mQ[3502] = aOmega*aCodonFreq[25];
	mQ[3518] = aK*aOmega*aCodonFreq[41];
	mQ[3522] = aOmega*aCodonFreq[45];
	mQ[3526] = aOmega*aCodonFreq[49];
	mQ[3530] = aK*aOmega*aCodonFreq[53];
	mQ[3534] = -aOmega*aCodonFreq[10]-aOmega*aCodonFreq[25]-aK*aOmega*aCodonFreq[41]-aOmega*aCodonFreq[45]-aOmega*aCodonFreq[49]-aK*aOmega*aCodonFreq[53]-aK*aCodonFreq[58]-aCodonFreq[59]-aCodonFreq[60];
	mQ[3535] = aK*aCodonFreq[58];
	mQ[3536] = aCodonFreq[59];
	mQ[3537] = aCodonFreq[60];
	mQ[3549] = aOmega*aCodonFreq[11];
	mQ[3564] = aOmega*aCodonFreq[26];
	mQ[3580] = aK*aOmega*aCodonFreq[42];
	mQ[3584] = aOmega*aCodonFreq[46];
	mQ[3588] = aOmega*aCodonFreq[50];
	mQ[3592] = aK*aOmega*aCodonFreq[54];
	mQ[3595] = aK*aCodonFreq[57];
	mQ[3596] = -aOmega*aCodonFreq[11]-aOmega*aCodonFreq[26]-aK*aOmega*aCodonFreq[42]-aOmega*aCodonFreq[46]-aOmega*aCodonFreq[50]-aK*aOmega*aCodonFreq[54]-aK*aCodonFreq[57]-aCodonFreq[59]-aCodonFreq[60];
	mQ[3597] = aCodonFreq[59];
	mQ[3598] = aCodonFreq[60];
	mQ[3626] = aOmega*aCodonFreq[27];
	mQ[3642] = aK*aOmega*aCodonFreq[43];
	mQ[3646] = aOmega*aCodonFreq[47];
	mQ[3650] = aOmega*aCodonFreq[51];
	mQ[3654] = aK*aOmega*aCodonFreq[55];
	mQ[3656] = aCodonFreq[57];
	mQ[3657] = aCodonFreq[58];
	mQ[3658] = -aOmega*aCodonFreq[27]-aK*aOmega*aCodonFreq[43]-aOmega*aCodonFreq[47]-aOmega*aCodonFreq[51]-aK*aOmega*aCodonFreq[55]-aCodonFreq[57]-aCodonFreq[58]-aK*aCodonFreq[60];
	mQ[3659] = aK*aCodonFreq[60];
	mQ[3672] = aOmega*aCodonFreq[12];
	mQ[3688] = aOmega*aCodonFreq[28];
	mQ[3704] = aK*aOmega*aCodonFreq[44];
	mQ[3708] = aOmega*aCodonFreq[48];
	mQ[3712] = aOmega*aCodonFreq[52];
	mQ[3716] = aK*aOmega*aCodonFreq[56];
	mQ[3717] = aCodonFreq[57];
	mQ[3718] = aCodonFreq[58];
	mQ[3719] = aK*aCodonFreq[59];
	mQ[3720] = -aOmega*aCodonFreq[12]-aOmega*aCodonFreq[28]-aK*aOmega*aCodonFreq[44]-aOmega*aCodonFreq[48]-aOmega*aCodonFreq[52]-aK*aOmega*aCodonFreq[56]-aCodonFreq[57]-aCodonFreq[58]-aK*aCodonFreq[59];

	// Compute the scale factor
	double scale_q =     aK*aCodonFreq[0]*aCodonFreq[1]
	             + aOmega*aCodonFreq[0]*aCodonFreq[2]
	             + aOmega*aCodonFreq[0]*aCodonFreq[3]
	             + aK*aOmega*aCodonFreq[0]*aCodonFreq[4]
	             + aOmega*aCodonFreq[0]*aCodonFreq[8]
	             + aOmega*aCodonFreq[0]*aCodonFreq[10]
	             + aK*aOmega*aCodonFreq[0]*aCodonFreq[13]
	             + aOmega*aCodonFreq[0]*aCodonFreq[29]
	             + aOmega*aCodonFreq[0]*aCodonFreq[45]
	             + aOmega*aCodonFreq[1]*aCodonFreq[2]
	             + aOmega*aCodonFreq[1]*aCodonFreq[3]
	             + aK*aOmega*aCodonFreq[1]*aCodonFreq[5]
	             + aOmega*aCodonFreq[1]*aCodonFreq[9]
	             + aOmega*aCodonFreq[1]*aCodonFreq[11]
	             + aK*aOmega*aCodonFreq[1]*aCodonFreq[14]
	             + aOmega*aCodonFreq[1]*aCodonFreq[30]
	             + aOmega*aCodonFreq[1]*aCodonFreq[46]
	             + aK*aCodonFreq[2]*aCodonFreq[3]
	             + aK*aOmega*aCodonFreq[2]*aCodonFreq[6]
	             + aK*aCodonFreq[2]*aCodonFreq[15]
	             + aOmega*aCodonFreq[2]*aCodonFreq[31]
	             + aOmega*aCodonFreq[2]*aCodonFreq[47]
	             + aK*aOmega*aCodonFreq[3]*aCodonFreq[7]
	             + aOmega*aCodonFreq[3]*aCodonFreq[12]
	             + aK*aCodonFreq[3]*aCodonFreq[16]
	             + aOmega*aCodonFreq[3]*aCodonFreq[32]
	             + aOmega*aCodonFreq[3]*aCodonFreq[48]
	             + aK*aCodonFreq[4]*aCodonFreq[5]
	             +   aCodonFreq[4]*aCodonFreq[6]
	             +   aCodonFreq[4]*aCodonFreq[7]
	             + aOmega*aCodonFreq[4]*aCodonFreq[8]
	             + aOmega*aCodonFreq[4]*aCodonFreq[10]
	             + aK*aOmega*aCodonFreq[4]*aCodonFreq[17]
	             + aOmega*aCodonFreq[4]*aCodonFreq[33]
	             + aOmega*aCodonFreq[4]*aCodonFreq[49]
	             +   aCodonFreq[5]*aCodonFreq[6]
	             +   aCodonFreq[5]*aCodonFreq[7]
	             + aOmega*aCodonFreq[5]*aCodonFreq[9]
	             + aOmega*aCodonFreq[5]*aCodonFreq[11]
	             + aK*aOmega*aCodonFreq[5]*aCodonFreq[18]
	             + aOmega*aCodonFreq[5]*aCodonFreq[34]
	             + aOmega*aCodonFreq[5]*aCodonFreq[50]
	             + aK*aCodonFreq[6]*aCodonFreq[7]
	             + aK*aOmega*aCodonFreq[6]*aCodonFreq[19]
	             + aOmega*aCodonFreq[6]*aCodonFreq[35]
	             + aOmega*aCodonFreq[6]*aCodonFreq[51]
	             + aOmega*aCodonFreq[7]*aCodonFreq[12]
	             + aK*aOmega*aCodonFreq[7]*aCodonFreq[20]
	             + aOmega*aCodonFreq[7]*aCodonFreq[36]
	             + aOmega*aCodonFreq[7]*aCodonFreq[52]
	             + aK*aCodonFreq[8]*aCodonFreq[9]
	             + aK*aOmega*aCodonFreq[8]*aCodonFreq[10]
	             + aK*aOmega*aCodonFreq[8]*aCodonFreq[21]
	             + aOmega*aCodonFreq[8]*aCodonFreq[37]
	             + aOmega*aCodonFreq[8]*aCodonFreq[53]
	             + aK*aOmega*aCodonFreq[9]*aCodonFreq[11]
	             + aK*aOmega*aCodonFreq[9]*aCodonFreq[22]
	             + aOmega*aCodonFreq[9]*aCodonFreq[38]
	             + aOmega*aCodonFreq[9]*aCodonFreq[54]
	             + aK*aCodonFreq[10]*aCodonFreq[11]
	             + aOmega*aCodonFreq[10]*aCodonFreq[12]
	             + aK*aOmega*aCodonFreq[10]*aCodonFreq[25]
	             + aOmega*aCodonFreq[10]*aCodonFreq[41]
	             + aOmega*aCodonFreq[10]*aCodonFreq[57]
	             + aOmega*aCodonFreq[11]*aCodonFreq[12]
	             + aK*aOmega*aCodonFreq[11]*aCodonFreq[26]
	             + aOmega*aCodonFreq[11]*aCodonFreq[42]
	             + aOmega*aCodonFreq[11]*aCodonFreq[58]
	             + aK*aOmega*aCodonFreq[12]*aCodonFreq[28]
	             + aOmega*aCodonFreq[12]*aCodonFreq[44]
	             + aOmega*aCodonFreq[12]*aCodonFreq[60]
	             + aK*aCodonFreq[13]*aCodonFreq[14]
	             +   aCodonFreq[13]*aCodonFreq[15]
	             +   aCodonFreq[13]*aCodonFreq[16]
	             + aK*aOmega*aCodonFreq[13]*aCodonFreq[17]
	             + aOmega*aCodonFreq[13]*aCodonFreq[21]
	             + aOmega*aCodonFreq[13]*aCodonFreq[25]
	             + aOmega*aCodonFreq[13]*aCodonFreq[29]
	             + aOmega*aCodonFreq[13]*aCodonFreq[45]
	             +   aCodonFreq[14]*aCodonFreq[15]
	             +   aCodonFreq[14]*aCodonFreq[16]
	             + aK*aOmega*aCodonFreq[14]*aCodonFreq[18]
	             + aOmega*aCodonFreq[14]*aCodonFreq[22]
	             + aOmega*aCodonFreq[14]*aCodonFreq[26]
	             + aOmega*aCodonFreq[14]*aCodonFreq[30]
	             + aOmega*aCodonFreq[14]*aCodonFreq[46]
	             + aK*aCodonFreq[15]*aCodonFreq[16]
	             + aK*aOmega*aCodonFreq[15]*aCodonFreq[19]
	             + aOmega*aCodonFreq[15]*aCodonFreq[23]
	             + aOmega*aCodonFreq[15]*aCodonFreq[27]
	             + aOmega*aCodonFreq[15]*aCodonFreq[31]
	             + aOmega*aCodonFreq[15]*aCodonFreq[47]
	             + aK*aOmega*aCodonFreq[16]*aCodonFreq[20]
	             + aOmega*aCodonFreq[16]*aCodonFreq[24]
	             + aOmega*aCodonFreq[16]*aCodonFreq[28]
	             + aOmega*aCodonFreq[16]*aCodonFreq[32]
	             + aOmega*aCodonFreq[16]*aCodonFreq[48]
	             + aK*aCodonFreq[17]*aCodonFreq[18]
	             +   aCodonFreq[17]*aCodonFreq[19]
	             +   aCodonFreq[17]*aCodonFreq[20]
	             + aOmega*aCodonFreq[17]*aCodonFreq[21]
	             + aOmega*aCodonFreq[17]*aCodonFreq[25]
	             + aOmega*aCodonFreq[17]*aCodonFreq[33]
	             + aOmega*aCodonFreq[17]*aCodonFreq[49]
	             +   aCodonFreq[18]*aCodonFreq[19]
	             +   aCodonFreq[18]*aCodonFreq[20]
	             + aOmega*aCodonFreq[18]*aCodonFreq[22]
	             + aOmega*aCodonFreq[18]*aCodonFreq[26]
	             + aOmega*aCodonFreq[18]*aCodonFreq[34]
	             + aOmega*aCodonFreq[18]*aCodonFreq[50]
	             + aK*aCodonFreq[19]*aCodonFreq[20]
	             + aOmega*aCodonFreq[19]*aCodonFreq[23]
	             + aOmega*aCodonFreq[19]*aCodonFreq[27]
	             + aOmega*aCodonFreq[19]*aCodonFreq[35]
	             + aOmega*aCodonFreq[19]*aCodonFreq[51]
	             + aOmega*aCodonFreq[20]*aCodonFreq[24]
	             + aOmega*aCodonFreq[20]*aCodonFreq[28]
	             + aOmega*aCodonFreq[20]*aCodonFreq[36]
	             + aOmega*aCodonFreq[20]*aCodonFreq[52]
	             + aK*aCodonFreq[21]*aCodonFreq[22]
	             + aOmega*aCodonFreq[21]*aCodonFreq[23]
	             + aOmega*aCodonFreq[21]*aCodonFreq[24]
	             + aK*aOmega*aCodonFreq[21]*aCodonFreq[25]
	             + aOmega*aCodonFreq[21]*aCodonFreq[37]
	             + aOmega*aCodonFreq[21]*aCodonFreq[53]
	             + aOmega*aCodonFreq[22]*aCodonFreq[23]
	             + aOmega*aCodonFreq[22]*aCodonFreq[24]
	             + aK*aOmega*aCodonFreq[22]*aCodonFreq[26]
	             + aOmega*aCodonFreq[22]*aCodonFreq[38]
	             + aOmega*aCodonFreq[22]*aCodonFreq[54]
	             + aK*aCodonFreq[23]*aCodonFreq[24]
	             + aK*aOmega*aCodonFreq[23]*aCodonFreq[27]
	             + aOmega*aCodonFreq[23]*aCodonFreq[39]
	             + aOmega*aCodonFreq[23]*aCodonFreq[55]
	             + aK*aOmega*aCodonFreq[24]*aCodonFreq[28]
	             + aOmega*aCodonFreq[24]*aCodonFreq[40]
	             + aOmega*aCodonFreq[24]*aCodonFreq[56]
	             + aK*aCodonFreq[25]*aCodonFreq[26]
	             +   aCodonFreq[25]*aCodonFreq[27]
	             +   aCodonFreq[25]*aCodonFreq[28]
	             + aOmega*aCodonFreq[25]*aCodonFreq[41]
	             + aOmega*aCodonFreq[25]*aCodonFreq[57]
	             +   aCodonFreq[26]*aCodonFreq[27]
	             +   aCodonFreq[26]*aCodonFreq[28]
	             + aOmega*aCodonFreq[26]*aCodonFreq[42]
	             + aOmega*aCodonFreq[26]*aCodonFreq[58]
	             + aK*aCodonFreq[27]*aCodonFreq[28]
	             +   aCodonFreq[27]*aCodonFreq[43]
	             + aOmega*aCodonFreq[27]*aCodonFreq[59]
	             +   aCodonFreq[28]*aCodonFreq[44]
	             + aOmega*aCodonFreq[28]*aCodonFreq[60]
	             + aK*aCodonFreq[29]*aCodonFreq[30]
	             +   aCodonFreq[29]*aCodonFreq[31]
	             + aOmega*aCodonFreq[29]*aCodonFreq[32]
	             + aK*aOmega*aCodonFreq[29]*aCodonFreq[33]
	             + aOmega*aCodonFreq[29]*aCodonFreq[37]
	             + aOmega*aCodonFreq[29]*aCodonFreq[41]
	             + aK*aOmega*aCodonFreq[29]*aCodonFreq[45]
	             +   aCodonFreq[30]*aCodonFreq[31]
	             + aOmega*aCodonFreq[30]*aCodonFreq[32]
	             + aK*aOmega*aCodonFreq[30]*aCodonFreq[34]
	             + aOmega*aCodonFreq[30]*aCodonFreq[38]
	             + aOmega*aCodonFreq[30]*aCodonFreq[42]
	             + aK*aOmega*aCodonFreq[30]*aCodonFreq[46]
	             + aK*aOmega*aCodonFreq[31]*aCodonFreq[32]
	             + aK*aOmega*aCodonFreq[31]*aCodonFreq[35]
	             + aOmega*aCodonFreq[31]*aCodonFreq[39]
	             + aOmega*aCodonFreq[31]*aCodonFreq[43]
	             + aK*aOmega*aCodonFreq[31]*aCodonFreq[47]
	             + aK*aOmega*aCodonFreq[32]*aCodonFreq[36]
	             + aOmega*aCodonFreq[32]*aCodonFreq[40]
	             + aOmega*aCodonFreq[32]*aCodonFreq[44]
	             + aK*aOmega*aCodonFreq[32]*aCodonFreq[48]
	             + aK*aCodonFreq[33]*aCodonFreq[34]
	             +   aCodonFreq[33]*aCodonFreq[35]
	             +   aCodonFreq[33]*aCodonFreq[36]
	             + aOmega*aCodonFreq[33]*aCodonFreq[37]
	             + aOmega*aCodonFreq[33]*aCodonFreq[41]
	             + aK*aOmega*aCodonFreq[33]*aCodonFreq[49]
	             +   aCodonFreq[34]*aCodonFreq[35]
	             +   aCodonFreq[34]*aCodonFreq[36]
	             + aOmega*aCodonFreq[34]*aCodonFreq[38]
	             + aOmega*aCodonFreq[34]*aCodonFreq[42]
	             + aK*aOmega*aCodonFreq[34]*aCodonFreq[50]
	             + aK*aCodonFreq[35]*aCodonFreq[36]
	             + aOmega*aCodonFreq[35]*aCodonFreq[39]
	             + aOmega*aCodonFreq[35]*aCodonFreq[43]
	             + aK*aOmega*aCodonFreq[35]*aCodonFreq[51]
	             + aOmega*aCodonFreq[36]*aCodonFreq[40]
	             + aOmega*aCodonFreq[36]*aCodonFreq[44]
	             + aK*aOmega*aCodonFreq[36]*aCodonFreq[52]
	             + aK*aCodonFreq[37]*aCodonFreq[38]
	             + aOmega*aCodonFreq[37]*aCodonFreq[39]
	             + aOmega*aCodonFreq[37]*aCodonFreq[40]
	             + aK*aOmega*aCodonFreq[37]*aCodonFreq[41]
	             + aK*aOmega*aCodonFreq[37]*aCodonFreq[53]
	             + aOmega*aCodonFreq[38]*aCodonFreq[39]
	             + aOmega*aCodonFreq[38]*aCodonFreq[40]
	             + aK*aOmega*aCodonFreq[38]*aCodonFreq[42]
	             + aK*aOmega*aCodonFreq[38]*aCodonFreq[54]
	             + aK*aCodonFreq[39]*aCodonFreq[40]
	             + aK*aOmega*aCodonFreq[39]*aCodonFreq[43]
	             + aK*aOmega*aCodonFreq[39]*aCodonFreq[55]
	             + aK*aOmega*aCodonFreq[40]*aCodonFreq[44]
	             + aK*aOmega*aCodonFreq[40]*aCodonFreq[56]
	             + aK*aCodonFreq[41]*aCodonFreq[42]
	             + aOmega*aCodonFreq[41]*aCodonFreq[43]
	             + aOmega*aCodonFreq[41]*aCodonFreq[44]
	             + aK*aOmega*aCodonFreq[41]*aCodonFreq[57]
	             + aOmega*aCodonFreq[42]*aCodonFreq[43]
	             + aOmega*aCodonFreq[42]*aCodonFreq[44]
	             + aK*aOmega*aCodonFreq[42]*aCodonFreq[58]
	             + aK*aCodonFreq[43]*aCodonFreq[44]
	             + aK*aOmega*aCodonFreq[43]*aCodonFreq[59]
	             + aK*aOmega*aCodonFreq[44]*aCodonFreq[60]
	             + aK*aCodonFreq[45]*aCodonFreq[46]
	             +   aCodonFreq[45]*aCodonFreq[47]
	             +   aCodonFreq[45]*aCodonFreq[48]
	             + aK*aOmega*aCodonFreq[45]*aCodonFreq[49]
	             + aOmega*aCodonFreq[45]*aCodonFreq[53]
	             + aOmega*aCodonFreq[45]*aCodonFreq[57]
	             +   aCodonFreq[46]*aCodonFreq[47]
	             +   aCodonFreq[46]*aCodonFreq[48]
	             + aK*aOmega*aCodonFreq[46]*aCodonFreq[50]
	             + aOmega*aCodonFreq[46]*aCodonFreq[54]
	             + aOmega*aCodonFreq[46]*aCodonFreq[58]
	             + aK*aCodonFreq[47]*aCodonFreq[48]
	             + aK*aOmega*aCodonFreq[47]*aCodonFreq[51]
	             + aOmega*aCodonFreq[47]*aCodonFreq[55]
	             + aOmega*aCodonFreq[47]*aCodonFreq[59]
	             + aK*aOmega*aCodonFreq[48]*aCodonFreq[52]
	             + aOmega*aCodonFreq[48]*aCodonFreq[56]
	             + aOmega*aCodonFreq[48]*aCodonFreq[60]
	             + aK*aCodonFreq[49]*aCodonFreq[50]
	             +   aCodonFreq[49]*aCodonFreq[51]
	             +   aCodonFreq[49]*aCodonFreq[52]
	             + aOmega*aCodonFreq[49]*aCodonFreq[53]
	             + aOmega*aCodonFreq[49]*aCodonFreq[57]
	             +   aCodonFreq[50]*aCodonFreq[51]
	             +   aCodonFreq[50]*aCodonFreq[52]
	             + aOmega*aCodonFreq[50]*aCodonFreq[54]
	             + aOmega*aCodonFreq[50]*aCodonFreq[58]
	             + aK*aCodonFreq[51]*aCodonFreq[52]
	             + aOmega*aCodonFreq[51]*aCodonFreq[55]
	             + aOmega*aCodonFreq[51]*aCodonFreq[59]
	             + aOmega*aCodonFreq[52]*aCodonFreq[56]
	             + aOmega*aCodonFreq[52]*aCodonFreq[60]
	             + aK*aCodonFreq[53]*aCodonFreq[54]
	             + aOmega*aCodonFreq[53]*aCodonFreq[55]
	             + aOmega*aCodonFreq[53]*aCodonFreq[56]
	             + aK*aOmega*aCodonFreq[53]*aCodonFreq[57]
	             + aOmega*aCodonFreq[54]*aCodonFreq[55]
	             + aOmega*aCodonFreq[54]*aCodonFreq[56]
	             + aK*aOmega*aCodonFreq[54]*aCodonFreq[58]
	             + aK*aCodonFreq[55]*aCodonFreq[56]
	             + aK*aOmega*aCodonFreq[55]*aCodonFreq[59]
	             + aK*aOmega*aCodonFreq[56]*aCodonFreq[60]
	             + aK*aCodonFreq[57]*aCodonFreq[58]
	             +   aCodonFreq[57]*aCodonFreq[59]
	             +   aCodonFreq[57]*aCodonFreq[60]
	             +   aCodonFreq[58]*aCodonFreq[59]
	             +   aCodonFreq[58]*aCodonFreq[60]
	             + aK*aCodonFreq[59]*aCodonFreq[60];

	return scale_q*2.;
}

double TransitionMatrix::fillQ(double aK, const double* aCodonFreq)
{
	mQ[   0] = -aK*aCodonFreq[1]-aCodonFreq[2]-aCodonFreq[3]-aK*aCodonFreq[4]-aCodonFreq[8]-aCodonFreq[10]-aK*aCodonFreq[13]-aCodonFreq[29]-aCodonFreq[45];
	mQ[   1] = aK*aCodonFreq[1];
	mQ[   2] = aCodonFreq[2];
	mQ[   3] = aCodonFreq[3];
	mQ[   4] = aK*aCodonFreq[4];
	mQ[   8] = aCodonFreq[8];
	mQ[  10] = aCodonFreq[10];
	mQ[  13] = aK*aCodonFreq[13];
	mQ[  29] = aCodonFreq[29];
	mQ[  45] = aCodonFreq[45];
	mQ[  61] = aK*aCodonFreq[0];
	mQ[  62] = -aK*aCodonFreq[0]-aCodonFreq[2]-aCodonFreq[3]-aK*aCodonFreq[5]-aCodonFreq[9]-aCodonFreq[11]-aK*aCodonFreq[14]-aCodonFreq[30]-aCodonFreq[46];
	mQ[  63] = aCodonFreq[2];
	mQ[  64] = aCodonFreq[3];
	mQ[  66] = aK*aCodonFreq[5];
	mQ[  70] = aCodonFreq[9];
	mQ[  72] = aCodonFreq[11];
	mQ[  75] = aK*aCodonFreq[14];
	mQ[  91] = aCodonFreq[30];
	mQ[ 107] = aCodonFreq[46];
	mQ[ 122] = aCodonFreq[0];
	mQ[ 123] = aCodonFreq[1];
	mQ[ 124] = -aCodonFreq[0]-aCodonFreq[1]-aK*aCodonFreq[3]-aK*aCodonFreq[6]-aK*aCodonFreq[15]-aCodonFreq[31]-aCodonFreq[47];
	mQ[ 125] = aK*aCodonFreq[3];
	mQ[ 128] = aK*aCodonFreq[6];
	mQ[ 137] = aK*aCodonFreq[15];
	mQ[ 153] = aCodonFreq[31];
	mQ[ 169] = aCodonFreq[47];
	mQ[ 183] = aCodonFreq[0];
	mQ[ 184] = aCodonFreq[1];
	mQ[ 185] = aK*aCodonFreq[2];
	mQ[ 186] = -aCodonFreq[0]-aCodonFreq[1]-aK*aCodonFreq[2]-aK*aCodonFreq[7]-aCodonFreq[12]-aK*aCodonFreq[16]-aCodonFreq[32]-aCodonFreq[48];
	mQ[ 190] = aK*aCodonFreq[7];
	mQ[ 195] = aCodonFreq[12];
	mQ[ 199] = aK*aCodonFreq[16];
	mQ[ 215] = aCodonFreq[32];
	mQ[ 231] = aCodonFreq[48];
	mQ[ 244] = aK*aCodonFreq[0];
	mQ[ 248] = -aK*aCodonFreq[0]-aK*aCodonFreq[5]-aCodonFreq[6]-aCodonFreq[7]-aCodonFreq[8]-aCodonFreq[10]-aK*aCodonFreq[17]-aCodonFreq[33]-aCodonFreq[49];
	mQ[ 249] = aK*aCodonFreq[5];
	mQ[ 250] = aCodonFreq[6];
	mQ[ 251] = aCodonFreq[7];
	mQ[ 252] = aCodonFreq[8];
	mQ[ 254] = aCodonFreq[10];
	mQ[ 261] = aK*aCodonFreq[17];
	mQ[ 277] = aCodonFreq[33];
	mQ[ 293] = aCodonFreq[49];
	mQ[ 306] = aK*aCodonFreq[1];
	mQ[ 309] = aK*aCodonFreq[4];
	mQ[ 310] = -aK*aCodonFreq[1]-aK*aCodonFreq[4]-aCodonFreq[6]-aCodonFreq[7]-aCodonFreq[9]-aCodonFreq[11]-aK*aCodonFreq[18]-aCodonFreq[34]-aCodonFreq[50];
	mQ[ 311] = aCodonFreq[6];
	mQ[ 312] = aCodonFreq[7];
	mQ[ 314] = aCodonFreq[9];
	mQ[ 316] = aCodonFreq[11];
	mQ[ 323] = aK*aCodonFreq[18];
	mQ[ 339] = aCodonFreq[34];
	mQ[ 355] = aCodonFreq[50];
	mQ[ 368] = aK*aCodonFreq[2];
	mQ[ 370] = aCodonFreq[4];
	mQ[ 371] = aCodonFreq[5];
	mQ[ 372] = -aK*aCodonFreq[2]-aCodonFreq[4]-aCodonFreq[5]-aK*aCodonFreq[7]-aK*aCodonFreq[19]-aCodonFreq[35]-aCodonFreq[51];
	mQ[ 373] = aK*aCodonFreq[7];
	mQ[ 385] = aK*aCodonFreq[19];
	mQ[ 401] = aCodonFreq[35];
	mQ[ 417] = aCodonFreq[51];
	mQ[ 430] = aK*aCodonFreq[3];
	mQ[ 431] = aCodonFreq[4];
	mQ[ 432] = aCodonFreq[5];
	mQ[ 433] = aK*aCodonFreq[6];
	mQ[ 434] = -aK*aCodonFreq[3]-aCodonFreq[4]-aCodonFreq[5]-aK*aCodonFreq[6]-aCodonFreq[12]-aK*aCodonFreq[20]-aCodonFreq[36]-aCodonFreq[52];
	mQ[ 439] = aCodonFreq[12];
	mQ[ 447] = aK*aCodonFreq[20];
	mQ[ 463] = aCodonFreq[36];
	mQ[ 479] = aCodonFreq[52];
	mQ[ 488] = aCodonFreq[0];
	mQ[ 492] = aCodonFreq[4];
	mQ[ 496] = -aCodonFreq[0]-aCodonFreq[4]-aK*aCodonFreq[9]-aK*aCodonFreq[10]-aK*aCodonFreq[21]-aCodonFreq[37]-aCodonFreq[53];
	mQ[ 497] = aK*aCodonFreq[9];
	mQ[ 498] = aK*aCodonFreq[10];
	mQ[ 509] = aK*aCodonFreq[21];
	mQ[ 525] = aCodonFreq[37];
	mQ[ 541] = aCodonFreq[53];
	mQ[ 550] = aCodonFreq[1];
	mQ[ 554] = aCodonFreq[5];
	mQ[ 557] = aK*aCodonFreq[8];
	mQ[ 558] = -aCodonFreq[1]-aCodonFreq[5]-aK*aCodonFreq[8]-aK*aCodonFreq[11]-aK*aCodonFreq[22]-aCodonFreq[38]-aCodonFreq[54];
	mQ[ 560] = aK*aCodonFreq[11];
	mQ[ 571] = aK*aCodonFreq[22];
	mQ[ 587] = aCodonFreq[38];
	mQ[ 603] = aCodonFreq[54];
	mQ[ 610] = aCodonFreq[0];
	mQ[ 614] = aCodonFreq[4];
	mQ[ 618] = aK*aCodonFreq[8];
	mQ[ 620] = -aCodonFreq[0]-aCodonFreq[4]-aK*aCodonFreq[8]-aK*aCodonFreq[11]-aCodonFreq[12]-aK*aCodonFreq[25]-aCodonFreq[41]-aCodonFreq[57];
	mQ[ 621] = aK*aCodonFreq[11];
	mQ[ 622] = aCodonFreq[12];
	mQ[ 635] = aK*aCodonFreq[25];
	mQ[ 651] = aCodonFreq[41];
	mQ[ 667] = aCodonFreq[57];
	mQ[ 672] = aCodonFreq[1];
	mQ[ 676] = aCodonFreq[5];
	mQ[ 680] = aK*aCodonFreq[9];
	mQ[ 681] = aK*aCodonFreq[10];
	mQ[ 682] = -aCodonFreq[1]-aCodonFreq[5]-aK*aCodonFreq[9]-aK*aCodonFreq[10]-aCodonFreq[12]-aK*aCodonFreq[26]-aCodonFreq[42]-aCodonFreq[58];
	mQ[ 683] = aCodonFreq[12];
	mQ[ 697] = aK*aCodonFreq[26];
	mQ[ 713] = aCodonFreq[42];
	mQ[ 729] = aCodonFreq[58];
	mQ[ 735] = aCodonFreq[3];
	mQ[ 739] = aCodonFreq[7];
	mQ[ 742] = aCodonFreq[10];
	mQ[ 743] = aCodonFreq[11];
	mQ[ 744] = -aCodonFreq[3]-aCodonFreq[7]-aCodonFreq[10]-aCodonFreq[11]-aK*aCodonFreq[28]-aCodonFreq[44]-aCodonFreq[60];
	mQ[ 760] = aK*aCodonFreq[28];
	mQ[ 776] = aCodonFreq[44];
	mQ[ 792] = aCodonFreq[60];
	mQ[ 793] = aK*aCodonFreq[0];
	mQ[ 806] = -aK*aCodonFreq[0]-aK*aCodonFreq[14]-aCodonFreq[15]-aCodonFreq[16]-aK*aCodonFreq[17]-aCodonFreq[21]-aCodonFreq[25]-aCodonFreq[29]-aCodonFreq[45];
	mQ[ 807] = aK*aCodonFreq[14];
	mQ[ 808] = aCodonFreq[15];
	mQ[ 809] = aCodonFreq[16];
	mQ[ 810] = aK*aCodonFreq[17];
	mQ[ 814] = aCodonFreq[21];
	mQ[ 818] = aCodonFreq[25];
	mQ[ 822] = aCodonFreq[29];
	mQ[ 838] = aCodonFreq[45];
	mQ[ 855] = aK*aCodonFreq[1];
	mQ[ 867] = aK*aCodonFreq[13];
	mQ[ 868] = -aK*aCodonFreq[1]-aK*aCodonFreq[13]-aCodonFreq[15]-aCodonFreq[16]-aK*aCodonFreq[18]-aCodonFreq[22]-aCodonFreq[26]-aCodonFreq[30]-aCodonFreq[46];
	mQ[ 869] = aCodonFreq[15];
	mQ[ 870] = aCodonFreq[16];
	mQ[ 872] = aK*aCodonFreq[18];
	mQ[ 876] = aCodonFreq[22];
	mQ[ 880] = aCodonFreq[26];
	mQ[ 884] = aCodonFreq[30];
	mQ[ 900] = aCodonFreq[46];
	mQ[ 917] = aK*aCodonFreq[2];
	mQ[ 928] = aCodonFreq[13];
	mQ[ 929] = aCodonFreq[14];
	mQ[ 930] = -aK*aCodonFreq[2]-aCodonFreq[13]-aCodonFreq[14]-aK*aCodonFreq[16]-aK*aCodonFreq[19]-aCodonFreq[23]-aCodonFreq[27]-aCodonFreq[31]-aCodonFreq[47];
	mQ[ 931] = aK*aCodonFreq[16];
	mQ[ 934] = aK*aCodonFreq[19];
	mQ[ 938] = aCodonFreq[23];
	mQ[ 942] = aCodonFreq[27];
	mQ[ 946] = aCodonFreq[31];
	mQ[ 962] = aCodonFreq[47];
	mQ[ 979] = aK*aCodonFreq[3];
	mQ[ 989] = aCodonFreq[13];
	mQ[ 990] = aCodonFreq[14];
	mQ[ 991] = aK*aCodonFreq[15];
	mQ[ 992] = -aK*aCodonFreq[3]-aCodonFreq[13]-aCodonFreq[14]-aK*aCodonFreq[15]-aK*aCodonFreq[20]-aCodonFreq[24]-aCodonFreq[28]-aCodonFreq[32]-aCodonFreq[48];
	mQ[ 996] = aK*aCodonFreq[20];
	mQ[1000] = aCodonFreq[24];
	mQ[1004] = aCodonFreq[28];
	mQ[1008] = aCodonFreq[32];
	mQ[1024] = aCodonFreq[48];
	mQ[1041] = aK*aCodonFreq[4];
	mQ[1050] = aK*aCodonFreq[13];
	mQ[1054] = -aK*aCodonFreq[4]-aK*aCodonFreq[13]-aK*aCodonFreq[18]-aCodonFreq[19]-aCodonFreq[20]-aCodonFreq[21]-aCodonFreq[25]-aCodonFreq[33]-aCodonFreq[49];
	mQ[1055] = aK*aCodonFreq[18];
	mQ[1056] = aCodonFreq[19];
	mQ[1057] = aCodonFreq[20];
	mQ[1058] = aCodonFreq[21];
	mQ[1062] = aCodonFreq[25];
	mQ[1070] = aCodonFreq[33];
	mQ[1086] = aCodonFreq[49];
	mQ[1103] = aK*aCodonFreq[5];
	mQ[1112] = aK*aCodonFreq[14];
	mQ[1115] = aK*aCodonFreq[17];
	mQ[1116] = -aK*aCodonFreq[5]-aK*aCodonFreq[14]-aK*aCodonFreq[17]-aCodonFreq[19]-aCodonFreq[20]-aCodonFreq[22]-aCodonFreq[26]-aCodonFreq[34]-aCodonFreq[50];
	mQ[1117] = aCodonFreq[19];
	mQ[1118] = aCodonFreq[20];
	mQ[1120] = aCodonFreq[22];
	mQ[1124] = aCodonFreq[26];
	mQ[1132] = aCodonFreq[34];
	mQ[1148] = aCodonFreq[50];
	mQ[1165] = aK*aCodonFreq[6];
	mQ[1174] = aK*aCodonFreq[15];
	mQ[1176] = aCodonFreq[17];
	mQ[1177] = aCodonFreq[18];
	mQ[1178] = -aK*aCodonFreq[6]-aK*aCodonFreq[15]-aCodonFreq[17]-aCodonFreq[18]-aK*aCodonFreq[20]-aCodonFreq[23]-aCodonFreq[27]-aCodonFreq[35]-aCodonFreq[51];
	mQ[1179] = aK*aCodonFreq[20];
	mQ[1182] = aCodonFreq[23];
	mQ[1186] = aCodonFreq[27];
	mQ[1194] = aCodonFreq[35];
	mQ[1210] = aCodonFreq[51];
	mQ[1227] = aK*aCodonFreq[7];
	mQ[1236] = aK*aCodonFreq[16];
	mQ[1237] = aCodonFreq[17];
	mQ[1238] = aCodonFreq[18];
	mQ[1239] = aK*aCodonFreq[19];
	mQ[1240] = -aK*aCodonFreq[7]-aK*aCodonFreq[16]-aCodonFreq[17]-aCodonFreq[18]-aK*aCodonFreq[19]-aCodonFreq[24]-aCodonFreq[28]-aCodonFreq[36]-aCodonFreq[52];
	mQ[1244] = aCodonFreq[24];
	mQ[1248] = aCodonFreq[28];
	mQ[1256] = aCodonFreq[36];
	mQ[1272] = aCodonFreq[52];
	mQ[1289] = aK*aCodonFreq[8];
	mQ[1294] = aCodonFreq[13];
	mQ[1298] = aCodonFreq[17];
	mQ[1302] = -aK*aCodonFreq[8]-aCodonFreq[13]-aCodonFreq[17]-aK*aCodonFreq[22]-aCodonFreq[23]-aCodonFreq[24]-aK*aCodonFreq[25]-aCodonFreq[37]-aCodonFreq[53];
	mQ[1303] = aK*aCodonFreq[22];
	mQ[1304] = aCodonFreq[23];
	mQ[1305] = aCodonFreq[24];
	mQ[1306] = aK*aCodonFreq[25];
	mQ[1318] = aCodonFreq[37];
	mQ[1334] = aCodonFreq[53];
	mQ[1351] = aK*aCodonFreq[9];
	mQ[1356] = aCodonFreq[14];
	mQ[1360] = aCodonFreq[18];
	mQ[1363] = aK*aCodonFreq[21];
	mQ[1364] = -aK*aCodonFreq[9]-aCodonFreq[14]-aCodonFreq[18]-aK*aCodonFreq[21]-aCodonFreq[23]-aCodonFreq[24]-aK*aCodonFreq[26]-aCodonFreq[38]-aCodonFreq[54];
	mQ[1365] = aCodonFreq[23];
	mQ[1366] = aCodonFreq[24];
	mQ[1368] = aK*aCodonFreq[26];
	mQ[1380] = aCodonFreq[38];
	mQ[1396] = aCodonFreq[54];
	mQ[1418] = aCodonFreq[15];
	mQ[1422] = aCodonFreq[19];
	mQ[1424] = aCodonFreq[21];
	mQ[1425] = aCodonFreq[22];
	mQ[1426] = -aCodonFreq[15]-aCodonFreq[19]-aCodonFreq[21]-aCodonFreq[22]-aK*aCodonFreq[24]-aK*aCodonFreq[27]-aCodonFreq[39]-aCodonFreq[55];
	mQ[1427] = aK*aCodonFreq[24];
	mQ[1430] = aK*aCodonFreq[27];
	mQ[1442] = aCodonFreq[39];
	mQ[1458] = aCodonFreq[55];
	mQ[1480] = aCodonFreq[16];
	mQ[1484] = aCodonFreq[20];
	mQ[1485] = aCodonFreq[21];
	mQ[1486] = aCodonFreq[22];
	mQ[1487] = aK*aCodonFreq[23];
	mQ[1488] = -aCodonFreq[16]-aCodonFreq[20]-aCodonFreq[21]-aCodonFreq[22]-aK*aCodonFreq[23]-aK*aCodonFreq[28]-aCodonFreq[40]-aCodonFreq[56];
	mQ[1492] = aK*aCodonFreq[28];
	mQ[1504] = aCodonFreq[40];
	mQ[1520] = aCodonFreq[56];
	mQ[1535] = aK*aCodonFreq[10];
	mQ[1538] = aCodonFreq[13];
	mQ[1542] = aCodonFreq[17];
	mQ[1546] = aK*aCodonFreq[21];
	mQ[1550] = -aK*aCodonFreq[10]-aCodonFreq[13]-aCodonFreq[17]-aK*aCodonFreq[21]-aK*aCodonFreq[26]-aCodonFreq[27]-aCodonFreq[28]-aCodonFreq[41]-aCodonFreq[57];
	mQ[1551] = aK*aCodonFreq[26];
	mQ[1552] = aCodonFreq[27];
	mQ[1553] = aCodonFreq[28];
	mQ[1566] = aCodonFreq[41];
	mQ[1582] = aCodonFreq[57];
	mQ[1597] = aK*aCodonFreq[11];
	mQ[1600] = aCodonFreq[14];
	mQ[1604] = aCodonFreq[18];
	mQ[1608] = aK*aCodonFreq[22];
	mQ[1611] = aK*aCodonFreq[25];
	mQ[1612] = -aK*aCodonFreq[11]-aCodonFreq[14]-aCodonFreq[18]-aK*aCodonFreq[22]-aK*aCodonFreq[25]-aCodonFreq[27]-aCodonFreq[28]-aCodonFreq[42]-aCodonFreq[58];
	mQ[1613] = aCodonFreq[27];
	mQ[1614] = aCodonFreq[28];
	mQ[1628] = aCodonFreq[42];
	mQ[1644] = aCodonFreq[58];
	mQ[1662] = aCodonFreq[15];
	mQ[1666] = aCodonFreq[19];
	mQ[1670] = aK*aCodonFreq[23];
	mQ[1672] = aCodonFreq[25];
	mQ[1673] = aCodonFreq[26];
	mQ[1674] = -aCodonFreq[15]-aCodonFreq[19]-aK*aCodonFreq[23]-aCodonFreq[25]-aCodonFreq[26]-aK*aCodonFreq[28]-aCodonFreq[43]-aCodonFreq[59];
	mQ[1675] = aK*aCodonFreq[28];
	mQ[1690] = aCodonFreq[43];
	mQ[1706] = aCodonFreq[59];
	mQ[1720] = aK*aCodonFreq[12];
	mQ[1724] = aCodonFreq[16];
	mQ[1728] = aCodonFreq[20];
	mQ[1732] = aK*aCodonFreq[24];
	mQ[1733] = aCodonFreq[25];
	mQ[1734] = aCodonFreq[26];
	mQ[1735] = aK*aCodonFreq[27];
	mQ[1736] = -aK*aCodonFreq[12]-aCodonFreq[16]-aCodonFreq[20]-aK*aCodonFreq[24]-aCodonFreq[25]-aCodonFreq[26]-aK*aCodonFreq[27]-aCodonFreq[44]-aCodonFreq[60];
	mQ[1752] = aCodonFreq[44];
	mQ[1768] = aCodonFreq[60];
	mQ[1769] = aCodonFreq[0];
	mQ[1782] = aCodonFreq[13];
	mQ[1798] = -aCodonFreq[0]-aCodonFreq[13]-aK*aCodonFreq[30]-aCodonFreq[31]-aCodonFreq[32]-aK*aCodonFreq[33]-aCodonFreq[37]-aCodonFreq[41]-aK*aCodonFreq[45];
	mQ[1799] = aK*aCodonFreq[30];
	mQ[1800] = aCodonFreq[31];
	mQ[1801] = aCodonFreq[32];
	mQ[1802] = aK*aCodonFreq[33];
	mQ[1806] = aCodonFreq[37];
	mQ[1810] = aCodonFreq[41];
	mQ[1814] = aK*aCodonFreq[45];
	mQ[1831] = aCodonFreq[1];
	mQ[1844] = aCodonFreq[14];
	mQ[1859] = aK*aCodonFreq[29];
	mQ[1860] = -aCodonFreq[1]-aCodonFreq[14]-aK*aCodonFreq[29]-aCodonFreq[31]-aCodonFreq[32]-aK*aCodonFreq[34]-aCodonFreq[38]-aCodonFreq[42]-aK*aCodonFreq[46];
	mQ[1861] = aCodonFreq[31];
	mQ[1862] = aCodonFreq[32];
	mQ[1864] = aK*aCodonFreq[34];
	mQ[1868] = aCodonFreq[38];
	mQ[1872] = aCodonFreq[42];
	mQ[1876] = aK*aCodonFreq[46];
	mQ[1893] = aCodonFreq[2];
	mQ[1906] = aCodonFreq[15];
	mQ[1920] = aCodonFreq[29];
	mQ[1921] = aCodonFreq[30];
	mQ[1922] = -aCodonFreq[2]-aCodonFreq[15]-aCodonFreq[29]-aCodonFreq[30]-aK*aCodonFreq[32]-aK*aCodonFreq[35]-aCodonFreq[39]-aCodonFreq[43]-aK*aCodonFreq[47];
	mQ[1923] = aK*aCodonFreq[32];
	mQ[1926] = aK*aCodonFreq[35];
	mQ[1930] = aCodonFreq[39];
	mQ[1934] = aCodonFreq[43];
	mQ[1938] = aK*aCodonFreq[47];
	mQ[1955] = aCodonFreq[3];
	mQ[1968] = aCodonFreq[16];
	mQ[1981] = aCodonFreq[29];
	mQ[1982] = aCodonFreq[30];
	mQ[1983] = aK*aCodonFreq[31];
	mQ[1984] = -aCodonFreq[3]-aCodonFreq[16]-aCodonFreq[29]-aCodonFreq[30]-aK*aCodonFreq[31]-aK*aCodonFreq[36]-aCodonFreq[40]-aCodonFreq[44]-aK*aCodonFreq[48];
	mQ[1988] = aK*aCodonFreq[36];
	mQ[1992] = aCodonFreq[40];
	mQ[1996] = aCodonFreq[44];
	mQ[2000] = aK*aCodonFreq[48];
	mQ[2017] = aCodonFreq[4];
	mQ[2030] = aCodonFreq[17];
	mQ[2042] = aK*aCodonFreq[29];
	mQ[2046] = -aCodonFreq[4]-aCodonFreq[17]-aK*aCodonFreq[29]-aK*aCodonFreq[34]-aCodonFreq[35]-aCodonFreq[36]-aCodonFreq[37]-aCodonFreq[41]-aK*aCodonFreq[49];
	mQ[2047] = aK*aCodonFreq[34];
	mQ[2048] = aCodonFreq[35];
	mQ[2049] = aCodonFreq[36];
	mQ[2050] = aCodonFreq[37];
	mQ[2054] = aCodonFreq[41];
	mQ[2062] = aK*aCodonFreq[49];
	mQ[2079] = aCodonFreq[5];
	mQ[2092] = aCodonFreq[18];
	mQ[2104] = aK*aCodonFreq[30];
	mQ[2107] = aK*aCodonFreq[33];
	mQ[2108] = -aCodonFreq[5]-aCodonFreq[18]-aK*aCodonFreq[30]-aK*aCodonFreq[33]-aCodonFreq[35]-aCodonFreq[36]-aCodonFreq[38]-aCodonFreq[42]-aK*aCodonFreq[50];
	mQ[2109] = aCodonFreq[35];
	mQ[2110] = aCodonFreq[36];
	mQ[2112] = aCodonFreq[38];
	mQ[2116] = aCodonFreq[42];
	mQ[2124] = aK*aCodonFreq[50];
	mQ[2141] = aCodonFreq[6];
	mQ[2154] = aCodonFreq[19];
	mQ[2166] = aK*aCodonFreq[31];
	mQ[2168] = aCodonFreq[33];
	mQ[2169] = aCodonFreq[34];
	mQ[2170] = -aCodonFreq[6]-aCodonFreq[19]-aK*aCodonFreq[31]-aCodonFreq[33]-aCodonFreq[34]-aK*aCodonFreq[36]-aCodonFreq[39]-aCodonFreq[43]-aK*aCodonFreq[51];
	mQ[2171] = aK*aCodonFreq[36];
	mQ[2174] = aCodonFreq[39];
	mQ[2178] = aCodonFreq[43];
	mQ[2186] = aK*aCodonFreq[51];
	mQ[2203] = aCodonFreq[7];
	mQ[2216] = aCodonFreq[20];
	mQ[2228] = aK*aCodonFreq[32];
	mQ[2229] = aCodonFreq[33];
	mQ[2230] = aCodonFreq[34];
	mQ[2231] = aK*aCodonFreq[35];
	mQ[2232] = -aCodonFreq[7]-aCodonFreq[20]-aK*aCodonFreq[32]-aCodonFreq[33]-aCodonFreq[34]-aK*aCodonFreq[35]-aCodonFreq[40]-aCodonFreq[44]-aK*aCodonFreq[52];
	mQ[2236] = aCodonFreq[40];
	mQ[2240] = aCodonFreq[44];
	mQ[2248] = aK*aCodonFreq[52];
	mQ[2265] = aCodonFreq[8];
	mQ[2278] = aCodonFreq[21];
	mQ[2286] = aCodonFreq[29];
	mQ[2290] = aCodonFreq[33];
	mQ[2294] = -aCodonFreq[8]-aCodonFreq[21]-aCodonFreq[29]-aCodonFreq[33]-aK*aCodonFreq[38]-aCodonFreq[39]-aCodonFreq[40]-aK*aCodonFreq[41]-aK*aCodonFreq[53];
	mQ[2295] = aK*aCodonFreq[38];
	mQ[2296] = aCodonFreq[39];
	mQ[2297] = aCodonFreq[40];
	mQ[2298] = aK*aCodonFreq[41];
	mQ[2310] = aK*aCodonFreq[53];
	mQ[2327] = aCodonFreq[9];
	mQ[2340] = aCodonFreq[22];
	mQ[2348] = aCodonFreq[30];
	mQ[2352] = aCodonFreq[34];
	mQ[2355] = aK*aCodonFreq[37];
	mQ[2356] = -aCodonFreq[9]-aCodonFreq[22]-aCodonFreq[30]-aCodonFreq[34]-aK*aCodonFreq[37]-aCodonFreq[39]-aCodonFreq[40]-aK*aCodonFreq[42]-aK*aCodonFreq[54];
	mQ[2357] = aCodonFreq[39];
	mQ[2358] = aCodonFreq[40];
	mQ[2360] = aK*aCodonFreq[42];
	mQ[2372] = aK*aCodonFreq[54];
	mQ[2402] = aCodonFreq[23];
	mQ[2410] = aCodonFreq[31];
	mQ[2414] = aCodonFreq[35];
	mQ[2416] = aCodonFreq[37];
	mQ[2417] = aCodonFreq[38];
	mQ[2418] = -aCodonFreq[23]-aCodonFreq[31]-aCodonFreq[35]-aCodonFreq[37]-aCodonFreq[38]-aK*aCodonFreq[40]-aK*aCodonFreq[43]-aK*aCodonFreq[55];
	mQ[2419] = aK*aCodonFreq[40];
	mQ[2422] = aK*aCodonFreq[43];
	mQ[2434] = aK*aCodonFreq[55];
	mQ[2464] = aCodonFreq[24];
	mQ[2472] = aCodonFreq[32];
	mQ[2476] = aCodonFreq[36];
	mQ[2477] = aCodonFreq[37];
	mQ[2478] = aCodonFreq[38];
	mQ[2479] = aK*aCodonFreq[39];
	mQ[2480] = -aCodonFreq[24]-aCodonFreq[32]-aCodonFreq[36]-aCodonFreq[37]-aCodonFreq[38]-aK*aCodonFreq[39]-aK*aCodonFreq[44]-aK*aCodonFreq[56];
	mQ[2484] = aK*aCodonFreq[44];
	mQ[2496] = aK*aCodonFreq[56];
	mQ[2511] = aCodonFreq[10];
	mQ[2526] = aCodonFreq[25];
	mQ[2530] = aCodonFreq[29];
	mQ[2534] = aCodonFreq[33];
	mQ[2538] = aK*aCodonFreq[37];
	mQ[2542] = -aCodonFreq[10]-aCodonFreq[25]-aCodonFreq[29]-aCodonFreq[33]-aK*aCodonFreq[37]-aK*aCodonFreq[42]-aCodonFreq[43]-aCodonFreq[44]-aK*aCodonFreq[57];
	mQ[2543] = aK*aCodonFreq[42];
	mQ[2544] = aCodonFreq[43];
	mQ[2545] = aCodonFreq[44];
	mQ[2558] = aK*aCodonFreq[57];
	mQ[2573] = aCodonFreq[11];
	mQ[2588] = aCodonFreq[26];
	mQ[2592] = aCodonFreq[30];
	mQ[2596] = aCodonFreq[34];
	mQ[2600] = aK*aCodonFreq[38];
	mQ[2603] = aK*aCodonFreq[41];
	mQ[2604] = -aCodonFreq[11]-aCodonFreq[26]-aCodonFreq[30]-aCodonFreq[34]-aK*aCodonFreq[38]-aK*aCodonFreq[41]-aCodonFreq[43]-aCodonFreq[44]-aK*aCodonFreq[58];
	mQ[2605] = aCodonFreq[43];
	mQ[2606] = aCodonFreq[44];
	mQ[2620] = aK*aCodonFreq[58];
	mQ[2650] = aCodonFreq[27];
	mQ[2654] = aCodonFreq[31];
	mQ[2658] = aCodonFreq[35];
	mQ[2662] = aK*aCodonFreq[39];
	mQ[2664] = aCodonFreq[41];
	mQ[2665] = aCodonFreq[42];
	mQ[2666] = -aCodonFreq[27]-aCodonFreq[31]-aCodonFreq[35]-aK*aCodonFreq[39]-aCodonFreq[41]-aCodonFreq[42]-aK*aCodonFreq[44]-aK*aCodonFreq[59];
	mQ[2667] = aK*aCodonFreq[44];
	mQ[2682] = aK*aCodonFreq[59];
	mQ[2696] = aCodonFreq[12];
	mQ[2712] = aCodonFreq[28];
	mQ[2716] = aCodonFreq[32];
	mQ[2720] = aCodonFreq[36];
	mQ[2724] = aK*aCodonFreq[40];
	mQ[2725] = aCodonFreq[41];
	mQ[2726] = aCodonFreq[42];
	mQ[2727] = aK*aCodonFreq[43];
	mQ[2728] = -aCodonFreq[12]-aCodonFreq[28]-aCodonFreq[32]-aCodonFreq[36]-aK*aCodonFreq[40]-aCodonFreq[41]-aCodonFreq[42]-aK*aCodonFreq[43]-aK*aCodonFreq[60];
	mQ[2744] = aK*aCodonFreq[60];
	mQ[2745] = aCodonFreq[0];
	mQ[2758] = aCodonFreq[13];
	mQ[2774] = aK*aCodonFreq[29];
	mQ[2790] = -aCodonFreq[0]-aCodonFreq[13]-aK*aCodonFreq[29]-aK*aCodonFreq[46]-aCodonFreq[47]-aCodonFreq[48]-aK*aCodonFreq[49]-aCodonFreq[53]-aCodonFreq[57];
	mQ[2791] = aK*aCodonFreq[46];
	mQ[2792] = aCodonFreq[47];
	mQ[2793] = aCodonFreq[48];
	mQ[2794] = aK*aCodonFreq[49];
	mQ[2798] = aCodonFreq[53];
	mQ[2802] = aCodonFreq[57];
	mQ[2807] = aCodonFreq[1];
	mQ[2820] = aCodonFreq[14];
	mQ[2836] = aK*aCodonFreq[30];
	mQ[2851] = aK*aCodonFreq[45];
	mQ[2852] = -aCodonFreq[1]-aCodonFreq[14]-aK*aCodonFreq[30]-aK*aCodonFreq[45]-aCodonFreq[47]-aCodonFreq[48]-aK*aCodonFreq[50]-aCodonFreq[54]-aCodonFreq[58];
	mQ[2853] = aCodonFreq[47];
	mQ[2854] = aCodonFreq[48];
	mQ[2856] = aK*aCodonFreq[50];
	mQ[2860] = aCodonFreq[54];
	mQ[2864] = aCodonFreq[58];
	mQ[2869] = aCodonFreq[2];
	mQ[2882] = aCodonFreq[15];
	mQ[2898] = aK*aCodonFreq[31];
	mQ[2912] = aCodonFreq[45];
	mQ[2913] = aCodonFreq[46];
	mQ[2914] = -aCodonFreq[2]-aCodonFreq[15]-aK*aCodonFreq[31]-aCodonFreq[45]-aCodonFreq[46]-aK*aCodonFreq[48]-aK*aCodonFreq[51]-aCodonFreq[55]-aCodonFreq[59];
	mQ[2915] = aK*aCodonFreq[48];
	mQ[2918] = aK*aCodonFreq[51];
	mQ[2922] = aCodonFreq[55];
	mQ[2926] = aCodonFreq[59];
	mQ[2931] = aCodonFreq[3];
	mQ[2944] = aCodonFreq[16];
	mQ[2960] = aK*aCodonFreq[32];
	mQ[2973] = aCodonFreq[45];
	mQ[2974] = aCodonFreq[46];
	mQ[2975] = aK*aCodonFreq[47];
	mQ[2976] = -aCodonFreq[3]-aCodonFreq[16]-aK*aCodonFreq[32]-aCodonFreq[45]-aCodonFreq[46]-aK*aCodonFreq[47]-aK*aCodonFreq[52]-aCodonFreq[56]-aCodonFreq[60];
	mQ[2980] = aK*aCodonFreq[52];
	mQ[2984] = aCodonFreq[56];
	mQ[2988] = aCodonFreq[60];
	mQ[2993] = aCodonFreq[4];
	mQ[3006] = aCodonFreq[17];
	mQ[3022] = aK*aCodonFreq[33];
	mQ[3034] = aK*aCodonFreq[45];
	mQ[3038] = -aCodonFreq[4]-aCodonFreq[17]-aK*aCodonFreq[33]-aK*aCodonFreq[45]-aK*aCodonFreq[50]-aCodonFreq[51]-aCodonFreq[52]-aCodonFreq[53]-aCodonFreq[57];
	mQ[3039] = aK*aCodonFreq[50];
	mQ[3040] = aCodonFreq[51];
	mQ[3041] = aCodonFreq[52];
	mQ[3042] = aCodonFreq[53];
	mQ[3046] = aCodonFreq[57];
	mQ[3055] = aCodonFreq[5];
	mQ[3068] = aCodonFreq[18];
	mQ[3084] = aK*aCodonFreq[34];
	mQ[3096] = aK*aCodonFreq[46];
	mQ[3099] = aK*aCodonFreq[49];
	mQ[3100] = -aCodonFreq[5]-aCodonFreq[18]-aK*aCodonFreq[34]-aK*aCodonFreq[46]-aK*aCodonFreq[49]-aCodonFreq[51]-aCodonFreq[52]-aCodonFreq[54]-aCodonFreq[58];
	mQ[3101] = aCodonFreq[51];
	mQ[3102] = aCodonFreq[52];
	mQ[3104] = aCodonFreq[54];
	mQ[3108] = aCodonFreq[58];
	mQ[3117] = aCodonFreq[6];
	mQ[3130] = aCodonFreq[19];
	mQ[3146] = aK*aCodonFreq[35];
	mQ[3158] = aK*aCodonFreq[47];
	mQ[3160] = aCodonFreq[49];
	mQ[3161] = aCodonFreq[50];
	mQ[3162] = -aCodonFreq[6]-aCodonFreq[19]-aK*aCodonFreq[35]-aK*aCodonFreq[47]-aCodonFreq[49]-aCodonFreq[50]-aK*aCodonFreq[52]-aCodonFreq[55]-aCodonFreq[59];
	mQ[3163] = aK*aCodonFreq[52];
	mQ[3166] = aCodonFreq[55];
	mQ[3170] = aCodonFreq[59];
	mQ[3179] = aCodonFreq[7];
	mQ[3192] = aCodonFreq[20];
	mQ[3208] = aK*aCodonFreq[36];
	mQ[3220] = aK*aCodonFreq[48];
	mQ[3221] = aCodonFreq[49];
	mQ[3222] = aCodonFreq[50];
	mQ[3223] = aK*aCodonFreq[51];
	mQ[3224] = -aCodonFreq[7]-aCodonFreq[20]-aK*aCodonFreq[36]-aK*aCodonFreq[48]-aCodonFreq[49]-aCodonFreq[50]-aK*aCodonFreq[51]-aCodonFreq[56]-aCodonFreq[60];
	mQ[3228] = aCodonFreq[56];
	mQ[3232] = aCodonFreq[60];
	mQ[3241] = aCodonFreq[8];
	mQ[3254] = aCodonFreq[21];
	mQ[3270] = aK*aCodonFreq[37];
	mQ[3278] = aCodonFreq[45];
	mQ[3282] = aCodonFreq[49];
	mQ[3286] = -aCodonFreq[8]-aCodonFreq[21]-aK*aCodonFreq[37]-aCodonFreq[45]-aCodonFreq[49]-aK*aCodonFreq[54]-aCodonFreq[55]-aCodonFreq[56]-aK*aCodonFreq[57];
	mQ[3287] = aK*aCodonFreq[54];
	mQ[3288] = aCodonFreq[55];
	mQ[3289] = aCodonFreq[56];
	mQ[3290] = aK*aCodonFreq[57];
	mQ[3303] = aCodonFreq[9];
	mQ[3316] = aCodonFreq[22];
	mQ[3332] = aK*aCodonFreq[38];
	mQ[3340] = aCodonFreq[46];
	mQ[3344] = aCodonFreq[50];
	mQ[3347] = aK*aCodonFreq[53];
	mQ[3348] = -aCodonFreq[9]-aCodonFreq[22]-aK*aCodonFreq[38]-aCodonFreq[46]-aCodonFreq[50]-aK*aCodonFreq[53]-aCodonFreq[55]-aCodonFreq[56]-aK*aCodonFreq[58];
	mQ[3349] = aCodonFreq[55];
	mQ[3350] = aCodonFreq[56];
	mQ[3352] = aK*aCodonFreq[58];
	mQ[3378] = aCodonFreq[23];
	mQ[3394] = aK*aCodonFreq[39];
	mQ[3402] = aCodonFreq[47];
	mQ[3406] = aCodonFreq[51];
	mQ[3408] = aCodonFreq[53];
	mQ[3409] = aCodonFreq[54];
	mQ[3410] = -aCodonFreq[23]-aK*aCodonFreq[39]-aCodonFreq[47]-aCodonFreq[51]-aCodonFreq[53]-aCodonFreq[54]-aK*aCodonFreq[56]-aK*aCodonFreq[59];
	mQ[3411] = aK*aCodonFreq[56];
	mQ[3414] = aK*aCodonFreq[59];
	mQ[3440] = aCodonFreq[24];
	mQ[3456] = aK*aCodonFreq[40];
	mQ[3464] = aCodonFreq[48];
	mQ[3468] = aCodonFreq[52];
	mQ[3469] = aCodonFreq[53];
	mQ[3470] = aCodonFreq[54];
	mQ[3471] = aK*aCodonFreq[55];
	mQ[3472] = -aCodonFreq[24]-aK*aCodonFreq[40]-aCodonFreq[48]-aCodonFreq[52]-aCodonFreq[53]-aCodonFreq[54]-aK*aCodonFreq[55]-aK*aCodonFreq[60];
	mQ[3476] = aK*aCodonFreq[60];
	mQ[3487] = aCodonFreq[10];
	mQ[3502] = aCodonFreq[25];
	mQ[3518] = aK*aCodonFreq[41];
	mQ[3522] = aCodonFreq[45];
	mQ[3526] = aCodonFreq[49];
	mQ[3530] = aK*aCodonFreq[53];
	mQ[3534] = -aCodonFreq[10]-aCodonFreq[25]-aK*aCodonFreq[41]-aCodonFreq[45]-aCodonFreq[49]-aK*aCodonFreq[53]-aK*aCodonFreq[58]-aCodonFreq[59]-aCodonFreq[60];
	mQ[3535] = aK*aCodonFreq[58];
	mQ[3536] = aCodonFreq[59];
	mQ[3537] = aCodonFreq[60];
	mQ[3549] = aCodonFreq[11];
	mQ[3564] = aCodonFreq[26];
	mQ[3580] = aK*aCodonFreq[42];
	mQ[3584] = aCodonFreq[46];
	mQ[3588] = aCodonFreq[50];
	mQ[3592] = aK*aCodonFreq[54];
	mQ[3595] = aK*aCodonFreq[57];
	mQ[3596] = -aCodonFreq[11]-aCodonFreq[26]-aK*aCodonFreq[42]-aCodonFreq[46]-aCodonFreq[50]-aK*aCodonFreq[54]-aK*aCodonFreq[57]-aCodonFreq[59]-aCodonFreq[60];
	mQ[3597] = aCodonFreq[59];
	mQ[3598] = aCodonFreq[60];
	mQ[3626] = aCodonFreq[27];
	mQ[3642] = aK*aCodonFreq[43];
	mQ[3646] = aCodonFreq[47];
	mQ[3650] = aCodonFreq[51];
	mQ[3654] = aK*aCodonFreq[55];
	mQ[3656] = aCodonFreq[57];
	mQ[3657] = aCodonFreq[58];
	mQ[3658] = -aCodonFreq[27]-aK*aCodonFreq[43]-aCodonFreq[47]-aCodonFreq[51]-aK*aCodonFreq[55]-aCodonFreq[57]-aCodonFreq[58]-aK*aCodonFreq[60];
	mQ[3659] = aK*aCodonFreq[60];
	mQ[3672] = aCodonFreq[12];
	mQ[3688] = aCodonFreq[28];
	mQ[3704] = aK*aCodonFreq[44];
	mQ[3708] = aCodonFreq[48];
	mQ[3712] = aCodonFreq[52];
	mQ[3716] = aK*aCodonFreq[56];
	mQ[3717] = aCodonFreq[57];
	mQ[3718] = aCodonFreq[58];
	mQ[3719] = aK*aCodonFreq[59];
	mQ[3720] = -aCodonFreq[12]-aCodonFreq[28]-aK*aCodonFreq[44]-aCodonFreq[48]-aCodonFreq[52]-aK*aCodonFreq[56]-aCodonFreq[57]-aCodonFreq[58]-aK*aCodonFreq[59];

	// Compute the scale factor
	double scale_q =     aK*aCodonFreq[0]*aCodonFreq[1]
	             +   aCodonFreq[0]*aCodonFreq[2]
	             +   aCodonFreq[0]*aCodonFreq[3]
	             + aK*aCodonFreq[0]*aCodonFreq[4]
	             +   aCodonFreq[0]*aCodonFreq[8]
	             +   aCodonFreq[0]*aCodonFreq[10]
	             + aK*aCodonFreq[0]*aCodonFreq[13]
	             +   aCodonFreq[0]*aCodonFreq[29]
	             +   aCodonFreq[0]*aCodonFreq[45]
	             +   aCodonFreq[1]*aCodonFreq[2]
	             +   aCodonFreq[1]*aCodonFreq[3]
	             + aK*aCodonFreq[1]*aCodonFreq[5]
	             +   aCodonFreq[1]*aCodonFreq[9]
	             +   aCodonFreq[1]*aCodonFreq[11]
	             + aK*aCodonFreq[1]*aCodonFreq[14]
	             +   aCodonFreq[1]*aCodonFreq[30]
	             +   aCodonFreq[1]*aCodonFreq[46]
	             + aK*aCodonFreq[2]*aCodonFreq[3]
	             + aK*aCodonFreq[2]*aCodonFreq[6]
	             + aK*aCodonFreq[2]*aCodonFreq[15]
	             +   aCodonFreq[2]*aCodonFreq[31]
	             +   aCodonFreq[2]*aCodonFreq[47]
	             + aK*aCodonFreq[3]*aCodonFreq[7]
	             +   aCodonFreq[3]*aCodonFreq[12]
	             + aK*aCodonFreq[3]*aCodonFreq[16]
	             +   aCodonFreq[3]*aCodonFreq[32]
	             +   aCodonFreq[3]*aCodonFreq[48]
	             + aK*aCodonFreq[4]*aCodonFreq[5]
	             +   aCodonFreq[4]*aCodonFreq[6]
	             +   aCodonFreq[4]*aCodonFreq[7]
	             +   aCodonFreq[4]*aCodonFreq[8]
	             +   aCodonFreq[4]*aCodonFreq[10]
	             + aK*aCodonFreq[4]*aCodonFreq[17]
	             +   aCodonFreq[4]*aCodonFreq[33]
	             +   aCodonFreq[4]*aCodonFreq[49]
	             +   aCodonFreq[5]*aCodonFreq[6]
	             +   aCodonFreq[5]*aCodonFreq[7]
	             +   aCodonFreq[5]*aCodonFreq[9]
	             +   aCodonFreq[5]*aCodonFreq[11]
	             + aK*aCodonFreq[5]*aCodonFreq[18]
	             +   aCodonFreq[5]*aCodonFreq[34]
	             +   aCodonFreq[5]*aCodonFreq[50]
	             + aK*aCodonFreq[6]*aCodonFreq[7]
	             + aK*aCodonFreq[6]*aCodonFreq[19]
	             +   aCodonFreq[6]*aCodonFreq[35]
	             +   aCodonFreq[6]*aCodonFreq[51]
	             +   aCodonFreq[7]*aCodonFreq[12]
	             + aK*aCodonFreq[7]*aCodonFreq[20]
	             +   aCodonFreq[7]*aCodonFreq[36]
	             +   aCodonFreq[7]*aCodonFreq[52]
	             + aK*aCodonFreq[8]*aCodonFreq[9]
	             + aK*aCodonFreq[8]*aCodonFreq[10]
	             + aK*aCodonFreq[8]*aCodonFreq[21]
	             +   aCodonFreq[8]*aCodonFreq[37]
	             +   aCodonFreq[8]*aCodonFreq[53]
	             + aK*aCodonFreq[9]*aCodonFreq[11]
	             + aK*aCodonFreq[9]*aCodonFreq[22]
	             +   aCodonFreq[9]*aCodonFreq[38]
	             +   aCodonFreq[9]*aCodonFreq[54]
	             + aK*aCodonFreq[10]*aCodonFreq[11]
	             +   aCodonFreq[10]*aCodonFreq[12]
	             + aK*aCodonFreq[10]*aCodonFreq[25]
	             +   aCodonFreq[10]*aCodonFreq[41]
	             +   aCodonFreq[10]*aCodonFreq[57]
	             +   aCodonFreq[11]*aCodonFreq[12]
	             + aK*aCodonFreq[11]*aCodonFreq[26]
	             +   aCodonFreq[11]*aCodonFreq[42]
	             +   aCodonFreq[11]*aCodonFreq[58]
	             + aK*aCodonFreq[12]*aCodonFreq[28]
	             +   aCodonFreq[12]*aCodonFreq[44]
	             +   aCodonFreq[12]*aCodonFreq[60]
	             + aK*aCodonFreq[13]*aCodonFreq[14]
	             +   aCodonFreq[13]*aCodonFreq[15]
	             +   aCodonFreq[13]*aCodonFreq[16]
	             + aK*aCodonFreq[13]*aCodonFreq[17]
	             +   aCodonFreq[13]*aCodonFreq[21]
	             +   aCodonFreq[13]*aCodonFreq[25]
	             +   aCodonFreq[13]*aCodonFreq[29]
	             +   aCodonFreq[13]*aCodonFreq[45]
	             +   aCodonFreq[14]*aCodonFreq[15]
	             +   aCodonFreq[14]*aCodonFreq[16]
	             + aK*aCodonFreq[14]*aCodonFreq[18]
	             +   aCodonFreq[14]*aCodonFreq[22]
	             +   aCodonFreq[14]*aCodonFreq[26]
	             +   aCodonFreq[14]*aCodonFreq[30]
	             +   aCodonFreq[14]*aCodonFreq[46]
	             + aK*aCodonFreq[15]*aCodonFreq[16]
	             + aK*aCodonFreq[15]*aCodonFreq[19]
	             +   aCodonFreq[15]*aCodonFreq[23]
	             +   aCodonFreq[15]*aCodonFreq[27]
	             +   aCodonFreq[15]*aCodonFreq[31]
	             +   aCodonFreq[15]*aCodonFreq[47]
	             + aK*aCodonFreq[16]*aCodonFreq[20]
	             +   aCodonFreq[16]*aCodonFreq[24]
	             +   aCodonFreq[16]*aCodonFreq[28]
	             +   aCodonFreq[16]*aCodonFreq[32]
	             +   aCodonFreq[16]*aCodonFreq[48]
	             + aK*aCodonFreq[17]*aCodonFreq[18]
	             +   aCodonFreq[17]*aCodonFreq[19]
	             +   aCodonFreq[17]*aCodonFreq[20]
	             +   aCodonFreq[17]*aCodonFreq[21]
	             +   aCodonFreq[17]*aCodonFreq[25]
	             +   aCodonFreq[17]*aCodonFreq[33]
	             +   aCodonFreq[17]*aCodonFreq[49]
	             +   aCodonFreq[18]*aCodonFreq[19]
	             +   aCodonFreq[18]*aCodonFreq[20]
	             +   aCodonFreq[18]*aCodonFreq[22]
	             +   aCodonFreq[18]*aCodonFreq[26]
	             +   aCodonFreq[18]*aCodonFreq[34]
	             +   aCodonFreq[18]*aCodonFreq[50]
	             + aK*aCodonFreq[19]*aCodonFreq[20]
	             +   aCodonFreq[19]*aCodonFreq[23]
	             +   aCodonFreq[19]*aCodonFreq[27]
	             +   aCodonFreq[19]*aCodonFreq[35]
	             +   aCodonFreq[19]*aCodonFreq[51]
	             +   aCodonFreq[20]*aCodonFreq[24]
	             +   aCodonFreq[20]*aCodonFreq[28]
	             +   aCodonFreq[20]*aCodonFreq[36]
	             +   aCodonFreq[20]*aCodonFreq[52]
	             + aK*aCodonFreq[21]*aCodonFreq[22]
	             +   aCodonFreq[21]*aCodonFreq[23]
	             +   aCodonFreq[21]*aCodonFreq[24]
	             + aK*aCodonFreq[21]*aCodonFreq[25]
	             +   aCodonFreq[21]*aCodonFreq[37]
	             +   aCodonFreq[21]*aCodonFreq[53]
	             +   aCodonFreq[22]*aCodonFreq[23]
	             +   aCodonFreq[22]*aCodonFreq[24]
	             + aK*aCodonFreq[22]*aCodonFreq[26]
	             +   aCodonFreq[22]*aCodonFreq[38]
	             +   aCodonFreq[22]*aCodonFreq[54]
	             + aK*aCodonFreq[23]*aCodonFreq[24]
	             + aK*aCodonFreq[23]*aCodonFreq[27]
	             +   aCodonFreq[23]*aCodonFreq[39]
	             +   aCodonFreq[23]*aCodonFreq[55]
	             + aK*aCodonFreq[24]*aCodonFreq[28]
	             +   aCodonFreq[24]*aCodonFreq[40]
	             +   aCodonFreq[24]*aCodonFreq[56]
	             + aK*aCodonFreq[25]*aCodonFreq[26]
	             +   aCodonFreq[25]*aCodonFreq[27]
	             +   aCodonFreq[25]*aCodonFreq[28]
	             +   aCodonFreq[25]*aCodonFreq[41]
	             +   aCodonFreq[25]*aCodonFreq[57]
	             +   aCodonFreq[26]*aCodonFreq[27]
	             +   aCodonFreq[26]*aCodonFreq[28]
	             +   aCodonFreq[26]*aCodonFreq[42]
	             +   aCodonFreq[26]*aCodonFreq[58]
	             + aK*aCodonFreq[27]*aCodonFreq[28]
	             +   aCodonFreq[27]*aCodonFreq[43]
	             +   aCodonFreq[27]*aCodonFreq[59]
	             +   aCodonFreq[28]*aCodonFreq[44]
	             +   aCodonFreq[28]*aCodonFreq[60]
	             + aK*aCodonFreq[29]*aCodonFreq[30]
	             +   aCodonFreq[29]*aCodonFreq[31]
	             +   aCodonFreq[29]*aCodonFreq[32]
	             + aK*aCodonFreq[29]*aCodonFreq[33]
	             +   aCodonFreq[29]*aCodonFreq[37]
	             +   aCodonFreq[29]*aCodonFreq[41]
	             + aK*aCodonFreq[29]*aCodonFreq[45]
	             +   aCodonFreq[30]*aCodonFreq[31]
	             +   aCodonFreq[30]*aCodonFreq[32]
	             + aK*aCodonFreq[30]*aCodonFreq[34]
	             +   aCodonFreq[30]*aCodonFreq[38]
	             +   aCodonFreq[30]*aCodonFreq[42]
	             + aK*aCodonFreq[30]*aCodonFreq[46]
	             + aK*aCodonFreq[31]*aCodonFreq[32]
	             + aK*aCodonFreq[31]*aCodonFreq[35]
	             +   aCodonFreq[31]*aCodonFreq[39]
	             +   aCodonFreq[31]*aCodonFreq[43]
	             + aK*aCodonFreq[31]*aCodonFreq[47]
	             + aK*aCodonFreq[32]*aCodonFreq[36]
	             +   aCodonFreq[32]*aCodonFreq[40]
	             +   aCodonFreq[32]*aCodonFreq[44]
	             + aK*aCodonFreq[32]*aCodonFreq[48]
	             + aK*aCodonFreq[33]*aCodonFreq[34]
	             +   aCodonFreq[33]*aCodonFreq[35]
	             +   aCodonFreq[33]*aCodonFreq[36]
	             +   aCodonFreq[33]*aCodonFreq[37]
	             +   aCodonFreq[33]*aCodonFreq[41]
	             + aK*aCodonFreq[33]*aCodonFreq[49]
	             +   aCodonFreq[34]*aCodonFreq[35]
	             +   aCodonFreq[34]*aCodonFreq[36]
	             +   aCodonFreq[34]*aCodonFreq[38]
	             +   aCodonFreq[34]*aCodonFreq[42]
	             + aK*aCodonFreq[34]*aCodonFreq[50]
	             + aK*aCodonFreq[35]*aCodonFreq[36]
	             +   aCodonFreq[35]*aCodonFreq[39]
	             +   aCodonFreq[35]*aCodonFreq[43]
	             + aK*aCodonFreq[35]*aCodonFreq[51]
	             +   aCodonFreq[36]*aCodonFreq[40]
	             +   aCodonFreq[36]*aCodonFreq[44]
	             + aK*aCodonFreq[36]*aCodonFreq[52]
	             + aK*aCodonFreq[37]*aCodonFreq[38]
	             +   aCodonFreq[37]*aCodonFreq[39]
	             +   aCodonFreq[37]*aCodonFreq[40]
	             + aK*aCodonFreq[37]*aCodonFreq[41]
	             + aK*aCodonFreq[37]*aCodonFreq[53]
	             +   aCodonFreq[38]*aCodonFreq[39]
	             +   aCodonFreq[38]*aCodonFreq[40]
	             + aK*aCodonFreq[38]*aCodonFreq[42]
	             + aK*aCodonFreq[38]*aCodonFreq[54]
	             + aK*aCodonFreq[39]*aCodonFreq[40]
	             + aK*aCodonFreq[39]*aCodonFreq[43]
	             + aK*aCodonFreq[39]*aCodonFreq[55]
	             + aK*aCodonFreq[40]*aCodonFreq[44]
	             + aK*aCodonFreq[40]*aCodonFreq[56]
	             + aK*aCodonFreq[41]*aCodonFreq[42]
	             +   aCodonFreq[41]*aCodonFreq[43]
	             +   aCodonFreq[41]*aCodonFreq[44]
	             + aK*aCodonFreq[41]*aCodonFreq[57]
	             +   aCodonFreq[42]*aCodonFreq[43]
	             +   aCodonFreq[42]*aCodonFreq[44]
	             + aK*aCodonFreq[42]*aCodonFreq[58]
	             + aK*aCodonFreq[43]*aCodonFreq[44]
	             + aK*aCodonFreq[43]*aCodonFreq[59]
	             + aK*aCodonFreq[44]*aCodonFreq[60]
	             + aK*aCodonFreq[45]*aCodonFreq[46]
	             +   aCodonFreq[45]*aCodonFreq[47]
	             +   aCodonFreq[45]*aCodonFreq[48]
	             + aK*aCodonFreq[45]*aCodonFreq[49]
	             +   aCodonFreq[45]*aCodonFreq[53]
	             +   aCodonFreq[45]*aCodonFreq[57]
	             +   aCodonFreq[46]*aCodonFreq[47]
	             +   aCodonFreq[46]*aCodonFreq[48]
	             + aK*aCodonFreq[46]*aCodonFreq[50]
	             +   aCodonFreq[46]*aCodonFreq[54]
	             +   aCodonFreq[46]*aCodonFreq[58]
	             + aK*aCodonFreq[47]*aCodonFreq[48]
	             + aK*aCodonFreq[47]*aCodonFreq[51]
	             +   aCodonFreq[47]*aCodonFreq[55]
	             +   aCodonFreq[47]*aCodonFreq[59]
	             + aK*aCodonFreq[48]*aCodonFreq[52]
	             +   aCodonFreq[48]*aCodonFreq[56]
	             +   aCodonFreq[48]*aCodonFreq[60]
	             + aK*aCodonFreq[49]*aCodonFreq[50]
	             +   aCodonFreq[49]*aCodonFreq[51]
	             +   aCodonFreq[49]*aCodonFreq[52]
	             +   aCodonFreq[49]*aCodonFreq[53]
	             +   aCodonFreq[49]*aCodonFreq[57]
	             +   aCodonFreq[50]*aCodonFreq[51]
	             +   aCodonFreq[50]*aCodonFreq[52]
	             +   aCodonFreq[50]*aCodonFreq[54]
	             +   aCodonFreq[50]*aCodonFreq[58]
	             + aK*aCodonFreq[51]*aCodonFreq[52]
	             +   aCodonFreq[51]*aCodonFreq[55]
	             +   aCodonFreq[51]*aCodonFreq[59]
	             +   aCodonFreq[52]*aCodonFreq[56]
	             +   aCodonFreq[52]*aCodonFreq[60]
	             + aK*aCodonFreq[53]*aCodonFreq[54]
	             +   aCodonFreq[53]*aCodonFreq[55]
	             +   aCodonFreq[53]*aCodonFreq[56]
	             + aK*aCodonFreq[53]*aCodonFreq[57]
	             +   aCodonFreq[54]*aCodonFreq[55]
	             +   aCodonFreq[54]*aCodonFreq[56]
	             + aK*aCodonFreq[54]*aCodonFreq[58]
	             + aK*aCodonFreq[55]*aCodonFreq[56]
	             + aK*aCodonFreq[55]*aCodonFreq[59]
	             + aK*aCodonFreq[56]*aCodonFreq[60]
	             + aK*aCodonFreq[57]*aCodonFreq[58]
	             +   aCodonFreq[57]*aCodonFreq[59]
	             +   aCodonFreq[57]*aCodonFreq[60]
	             +   aCodonFreq[58]*aCodonFreq[59]
	             +   aCodonFreq[58]*aCodonFreq[60]
	             + aK*aCodonFreq[59]*aCodonFreq[60];

	return scale_q*2.;
}

