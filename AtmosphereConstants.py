# -*- coding: utf-8 -*-
"""
Created on Mon May  9 09:39:54 2016

This file contains the constans that are used inside the atmospheric modeling
code in Atmosphere.py (defining the Atmosphere class). This file does not need
any editing if one is not attempting to edit or add a new model. When adding
a new model support it in the Atmosphere's __init__(...) function by defining
a new model name and link the class' atmospheric interpolation constants to the
ones defined in here.

@author: MaxHenger
"""

class AtmospherePreliminary:
    def __init__(self):
        self.tkcDeepTemp = [ [ 1000.0, 1000.0, 1000.0, 1000.0, 3000.0, 
                4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 
                9000.0, 10000.0, 11000.0, 12000.0, 13000.0, 
                14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 
                19000.0, 20000.0, 21000.0, 22000.0, 23000.0, 
                24000.0, 25000.0, 26000.0, 27000.0, 28000.0, 
                29000.0, 30000.0, 32000.0, 32000.0, 32000.0, 
                32000.0 ],
                [ 727.7, 722.878179663, 714.993640674, 704.584538989, 696.87047879, 
                688.733545849, 680.995337814, 673.885102894, 665.664250609, 658.257894668, 
                650.504170719, 643.325422457, 635.394139452, 628.098019736, 620.813781605, 
                613.446853845, 605.198803014, 596.957934097, 589.569460597, 580.564223517, 
                572.373645337, 564.341195136, 556.06157412, 547.412508383, 539.288392349, 
                530.633922221, 522.375918767, 513.662402711, 505.774470388, 493.861464186, 
                485.469267907, 479.9, 0.0, 0.0, 0.0, 
                0.0 ],
                3 ]
        self.tkcDeepPres = [ [ 1000.0, 1000.0, 1000.0, 1000.0, 3000.0, 
                4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 
                9000.0, 10000.0, 11000.0, 12000.0, 13000.0, 
                14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 
                19000.0, 20000.0, 21000.0, 22000.0, 23000.0, 
                24000.0, 25000.0, 26000.0, 27000.0, 28000.0, 
                29000.0, 30000.0, 32000.0, 32000.0, 32000.0, 
                32000.0 ],
                [ 8645000.0, 8278118.51302, 7756762.97395, 7115688.8724, 6660750.06182, 
                6231310.8803, 5824006.41697, 5440663.45181, 5077339.7758, 4735977.44501, 
                4412750.44416, 4109020.77834, 3823166.44249, 3554313.45169, 3301579.75076, 
                3063367.54527, 2840950.06817, 2630832.18204, 2433721.20368, 2250283.00325, 
                2077146.78332, 1915129.86347, 1764333.76281, 1623535.08528, 1491525.89608, 
                1368361.33041, 1255028.78228, 1147523.54046, 1048877.05588, 926331.962748, 
                844609.018626, 794000.0, 0.0, 0.0, 0.0, 
                0.0 ],
                3 ]
        self.tkcDeepDens = [ [ 1000.0, 1000.0, 1000.0, 1000.0, 3000.0, 
                4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 
                9000.0, 10000.0, 11000.0, 12000.0, 13000.0, 
                14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 
                19000.0, 20000.0, 21000.0, 22000.0, 23000.0, 
                24000.0, 25000.0, 26000.0, 27000.0, 28000.0, 
                29000.0, 30000.0, 32000.0, 32000.0, 32000.0, 
                32000.0 ],
                [ 61.56, 59.4460524481, 56.3928951039, 52.6064906776, 49.8477432925, 
                47.2225361525, 44.7021120976, 42.229015457, 39.9418260743, 37.7036802457, 
                35.5634529428, 33.5225079832, 31.5865151245, 29.7314315189, 27.9277588, 
                26.257533281, 24.6621080759, 23.1740344155, 21.7217542622, 20.3789485356, 
                19.1024515952, 17.8712450836, 16.6925680703, 15.6184826351, 14.5535013894, 
                13.5875118072, 12.6364513818, 11.7666826656, 10.9168179558, 9.88912136277, 
                9.15393931861, 8.704, 0.0, 0.0, 0.0, 
                0.0 ],
                3 ]
        self.tkcUpperTemp = [ [ 33000.0, 33000.0, 33000.0, 33000.0, 35000.0, 
                36000.0, 37000.0, 38000.0, 39000.0, 40000.0, 
                41000.0, 42000.0, 43000.0, 44000.0, 45000.0, 
                46000.0, 47000.0, 48000.0, 49000.0, 50000.0, 
                51000.0, 52000.0, 53000.0, 54000.0, 55000.0, 
                56000.0, 57000.0, 58000.0, 59000.0, 60000.0, 
                62000.0, 64000.0, 66000.0, 68000.0, 70000.0, 
                72000.0, 74000.0, 76000.0, 78000.0, 80000.0, 
                82000.0, 84000.0, 86000.0, 88000.0, 90000.0, 
                92000.0, 94000.0, 96000.0, 100000.0, 100000.0, 
                100000.0, 100000.0 ],
                [ 30.0, 30.0, 30.0, 30.0, 60.0, 
                85.0, 85.0, 85.0, 85.0 ],
                [ 471.7, 471.874193548, 471.235483871, 472.353225806, 470.2, 
                467.243912054, 466.274060683, 465.851849774, 466.931203769, 464.528467772, 
                455.512175892, 458.427900139, 456.763020883, 458.855543478, 455.543064457, 
                451.565069495, 445.945587068, 446.416284088, 447.687707284, 443.718736648, 
                427.214322725, 441.716781469, 436.907989329, 438.715335537, 434.738070902, 
                478.977639604, 417.500835442, 436.582296232, 429.543423686, 426.128979743, 
                251.875118861, 472.687618697, 397.842180582, 426.317779755, 416.746010126, 
                464.121884954, 402.520302671, 420.758013698, 411.856783458, 407.686979754, 
                397.237341325, 411.535686747, 405.713721616, 402.819244118, 399.306070859, 
                406.928749746, 400.237595503, 403.452046074, 393.532548315, 390.28873681, 
                396.047659692, 397.023608662, 391.418954302, 386.111852944, 381.938981901, 
                391.480611487, 388.40861501, 388.397082955, 375.530792598, 372.755335585, 
                385.229894361, 381.092899039, 382.890133231, 366.747772363, 363.23967576, 
                379.999811069, 376.454627542, 376.849480894, 356.669874222, 354.685961377, 
                372.970861364, 371.497623052, 368.087857172, 350.056601715, 344.416478733, 
                366.716743475, 367.005202832, 363.564896871, 339.531317485, 337.648123691, 
                358.562164734, 359.018339815, 357.354490829, 331.719920462, 327.591026502, 
                350.634597588, 350.615631457, 351.932623683, 323.458893142, 319.387770303, 
                341.899444915, 340.683005326, 344.544691858, 313.709023099, 309.457892288, 
                333.767622751, 333.58266982, 335.774415338, 306.865946362, 301.580660544, 
                322.830064081, 321.591476684, 323.610550018, 296.557657403, 290.219465536, 
                312.912120925, 314.923681507, 311.724029753, 287.241058433, 280.54147731, 
                302.32145222, 302.34476503, 302.010750323, 276.372015673, 270.014625222, 
                291.602070195, 292.502419665, 288.585872179, 268.745789269, 258.200021801, 
                282.070266999, 278.473943409, 278.970062037, 251.784612196, 247.585287573, 
                275.116861807, 270.666322829, 256.228503328, 262.879919654, 238.258827906, 
                268.662285772, 241.051733016, 301.273344005, 172.300727108, 238.179400801, 
                260.357898218, 351.008093225, -12.514982475, 585.380240133, 243.838291586, 
                255.46748718, 225.507205574, 305.531563108, 101.830475666, 252.058371725, 
                243.794070856, 267.572584368, 177.10556604, 272.349194632, 244.172237543, 
                241.756229396, 242.665682762, 217.377570581, 235.759842578, 240.652678105, 
                235.181011561, 236.100168457, 218.956194645, 249.480610681, 245.217050039, 
                229.919724362, 229.089772444, 222.281306752, 250.919865235, 243.679121741, 
                223.940090993, 224.406548217, 221.243094476, 251.792258128, 241.066462998, 
                218.919911666, 219.957583074, 221.916852978, 245.903575369, 241.055026267, 
                211.980262343, 214.936022714, 217.328418344, 242.152938602, 235.913431933, 
                205.759038961, 207.91251962, 218.331624183, 224.021229364, 225.291246001, 
                196.783581811, 200.997124612, 204.656482772, 218.749240718, 212.921584064, 
                189.706633795, 195.597046449, 196.580939351, 207.260660812, 206.422417743, 
                183.78988301, 188.427592817, 190.318684556, 196.900947575, 198.588744963, 
                177.933834164, 181.930001639, 180.944537477, 187.503649247, 190.822602403, 
                173.674780333, 176.814981272, 176.80295048, 182.549688411, 186.520845424, 
                168.967044504, 172.505557144, 170.255703613, 178.72232829, 180.894015899, 
                166.857041651, 171.032467573, 169.888428617, 175.944510973, 177.703090979, 
                166.804788893, 169.767153211, 168.517033532, 172.142638571, 173.293620183, 
                169.896807405, 170.403916975, 170.478260003, 168.446584242, 169.065364322, 
                172.901596298, 172.269278072, 168.533405841, 165.64453942, 166.592317839, 
                175.4, 173.421075269, 169.066021505, 165.147550777, 166.1 ],
                3,
                3 ]
        self.tkcUpperPres = [ [ 33000.0, 33000.0, 33000.0, 33000.0, 35000.0, 
                36000.0, 37000.0, 38000.0, 39000.0, 40000.0, 
                41000.0, 42000.0, 43000.0, 44000.0, 45000.0, 
                46000.0, 47000.0, 48000.0, 49000.0, 50000.0, 
                51000.0, 52000.0, 53000.0, 54000.0, 55000.0, 
                56000.0, 57000.0, 58000.0, 59000.0, 60000.0, 
                62000.0, 64000.0, 66000.0, 68000.0, 70000.0, 
                72000.0, 74000.0, 76000.0, 78000.0, 80000.0, 
                82000.0, 84000.0, 86000.0, 88000.0, 90000.0, 
                92000.0, 94000.0, 96000.0, 100000.0, 100000.0, 
                100000.0, 100000.0 ],
                [ 30.0, 30.0, 30.0, 30.0, 60.0, 
                85.0, 85.0, 85.0, 85.0 ],
                [ 721100.0, 721111.612903, 721069.032258, 721143.548387, 721000.0, 
                674235.876176, 653258.652711, 721245.82548, 662380.346627, 674026.67789, 
                609928.247648, 652051.189201, 515672.918933, 634082.735659, 609496.64422, 
                533807.628528, 470385.456342, 674505.666404, 497846.256926, 532796.700337, 
                481471.393021, 720109.87026, -49834.3090508, 618398.883203, 480091.540405, 
                433506.799389, 369276.998102, 574859.741842, 396757.134992, 431637.138045, 
                389701.409422, 406138.266366, 350612.33093, 399274.727367, 387159.907416, 
                349487.562922, 344546.710629, 358586.203256, 346045.747654, 346523.23229, 
                312948.338891, 313563.923374, 308572.103358, 312659.486316, 309347.163425, 
                279719.081514, 278411.789421, 279154.200518, 277697.310665, 275488.114009, 
                249375.335052, 249320.531845, 247193.460163, 247455.930523, 244700.380537, 
                222179.578276, 221091.889652, 220843.141627, 219142.408105, 216910.363843, 
                197506.351843, 196750.619226, 196064.080858, 193808.487237, 191458.16409, 
                175195.01435, 174350.794734, 173513.438169, 170776.331118, 168656.979795, 
                155313.590756, 154113.298611, 153869.908402, 150026.690081, 147913.916729, 
                137150.622627, 136327.623726, 135822.627148, 131432.679166, 129487.353287, 
                121083.918737, 120265.238744, 119268.830317, 114782.019779, 112936.670122, 
                106313.702425, 105743.034202, 105351.083841, 99932.7901049, 97965.9662262, 
                93261.2715625, 92568.5599304, 92177.6730271, 86536.8556436, 84759.4649736, 
                81461.2113249, 80882.7260759, 80711.5573834, 74703.1206541, 72996.1738793, 
                70913.8831378, 70367.2454433, 70244.8716327, 64236.9341414, 62575.839509, 
                61423.2561238, 60979.2598927, 60866.3754406, 54951.9384793, 53320.4680846, 
                52993.092367, 52608.1020827, 52476.594347, 46740.6524432, 45202.2881526, 
                45444.3744082, 45109.7511312, 44973.4622252, 39526.8854396, 38030.3793051, 
                38769.4100002, 38498.1837151, 38238.7825588, 33216.6445081, 31756.1946271, 
                32937.9855911, 32543.5785247, 32375.2354967, 27656.5006858, 26344.8421864, 
                27838.6476355, 27360.4054118, 27182.5335191, 22941.2954009, 21744.4366273, 
                22010.3492776, 21488.3291494, 21099.8240201, 17558.4004157, 16775.0695299, 
                16263.6183242, 15905.9138218, 14964.5076218, 12496.8619061, 12252.1910942, 
                11325.8312349, 11197.4679818, 9719.25425599, 8472.99121118, 8491.59040598, 
                7793.05673612, 7783.67231549, 6287.25384884, 5854.36321328, 5837.44728183, 
                5321.94182062, 5314.08791753, 4153.0765852, 4054.30862385, 4020.6204667, 
                3601.17598142, 3590.51794986, 2762.32798241, 2824.99010492, 2776.07085138, 
                2413.35425372, 2403.40157333, 1840.78137764, 1972.08579516, 1909.09612779, 
                1601.4070037, 1597.772531, 1226.48844252, 1368.83517321, 1311.54463746, 
                1051.01773147, 1051.85668977, 814.002486693, 943.602185824, 896.725322357, 
                680.522070424, 685.860064745, 537.64333114, 638.356441914, 604.75407311, 
                433.493986834, 440.741760925, 349.887629608, 423.539430032, 399.458385202, 
                271.101982239, 279.0954722, 222.679268709, 276.273293155, 259.212386081, 
                166.898084211, 174.086672856, 139.301102008, 176.328329248, 165.692070473, 
                101.105680917, 106.675255732, 85.4698716465, 109.837045767, 104.019332027, 
                60.4791921219, 64.4819816364, 51.4669382879, 67.384778005, 64.4306014184, 
                35.6975505957, 38.4566886902, 30.5427192878, 40.8739138971, 39.4182622994, 
                20.8906054952, 22.8169410223, 17.8937114427, 24.5853011735, 23.8363493841, 
                12.2000274235, 13.4674181886, 10.5241123608, 14.6549530934, 14.2763401642, 
                5.49237060655, 6.11984868671, 4.75537527201, 6.65978658103, 6.51893989051, 
                3.60906469672, 3.9946890975, 3.09495429948, 4.26782609419, 4.18128005475, 
                2.66, 2.93835698925, 2.25993691756, 3.10657945042, 3.048 ],
                3,
                3 ]
        self.tkcUpperDens = [ [ 33000.0, 33000.0, 33000.0, 33000.0, 35000.0, 
                36000.0, 37000.0, 38000.0, 39000.0, 40000.0, 
                41000.0, 42000.0, 43000.0, 44000.0, 45000.0, 
                46000.0, 47000.0, 48000.0, 49000.0, 50000.0, 
                51000.0, 52000.0, 53000.0, 54000.0, 55000.0, 
                56000.0, 57000.0, 58000.0, 59000.0, 60000.0, 
                62000.0, 64000.0, 66000.0, 68000.0, 70000.0, 
                72000.0, 74000.0, 76000.0, 78000.0, 80000.0, 
                82000.0, 84000.0, 86000.0, 88000.0, 90000.0, 
                92000.0, 94000.0, 96000.0, 100000.0, 100000.0, 
                100000.0, 100000.0 ],
                [ 30.0, 30.0, 30.0, 30.0, 60.0, 
                85.0, 85.0, 85.0, 85.0 ],
                [ 8.041, 8.03809677419, 8.04874193548, 8.03011290323, 8.066, 
                7.61493132048, 7.60907064209, 7.62216464477, 7.60043904597, 7.63866156959, 
                7.01413735904, 7.0120318341, 7.02157017283, 7.00446091642, 7.04367686081, 
                6.2644606281, 6.28634275064, 6.29570988413, 6.26785248239, 6.31831804212, 
                5.75943931128, 5.78521943419, 5.74485536283, 5.80739237362, 5.8083906754, 
                5.26978212677, 5.26388919003, 5.36424285809, 5.18207264678, 5.3261192563, 
                4.81743218162, 4.82157219278, 4.83324417255, 4.88745682422, 4.88513229941, 
                4.39848914674, 4.41503494208, 4.42121270979, 4.43098177678, 4.46735154607, 
                4.01261123143, 4.018636426, 4.02130928937, 4.07142252026, 4.07146151632, 
                3.64106592753, 3.65409677329, 3.6237436811, 3.71272957585, 3.70880238866, 
                3.29912505845, 3.30404744858, 3.30152673892, 3.36984197203, 3.36532892904, 
                2.98043383867, 2.9943585937, 2.98476226646, 3.06298855753, 3.05588189517, 
                2.68913958688, 2.70588591858, 2.68311021676, 2.77633283011, 2.76514349028, 
                2.42100781381, 2.43332353845, 2.41586138263, 2.51099911844, 2.50354414371, 
                2.18282915789, 2.17736186311, 2.18733242477, 2.25481406529, 2.24467993488, 
                1.96367555462, 1.94971933171, 1.96216805808, 2.02680555229, 2.02373611678, 
                1.76446862363, 1.75198016491, 1.73541039668, 1.82450852843, 1.76837559801, 
                1.59244995086, 1.5735729119, 1.60095594661, 1.56904205441, 1.73276149117, 
                1.42973157294, 1.42851528428, 1.3913335588, 1.45999708905, 1.40457843729, 
                1.28062375738, 1.26805627358, 1.26586895796, 1.27291941016, 1.27692475967, 
                1.15177339756, 1.15270478271, 1.13467017925, 1.13630018068, 1.13172252404, 
                1.0302826524, 1.01573104721, 1.01816645405, 1.01276517176, 0.988185144158, 
                0.919095992842, 0.906074254268, 0.930122069045, 0.862026229033, 0.915536899326, 
                0.817533376233, 0.833682903464, 0.752616022453, 0.86489048558, 0.63266725854, 
                0.720570502226, 0.632827680265, 0.96555771211, 0.347170968433, 1.19719406651, 
                0.627384614864, 0.972861859346, -0.25006149455, 1.8341837052, -1.37684352459, 
                0.543291038318, -0.681326730552, 3.8558259005, -4.12709933764, 7.77818003186, 
                0.443237963046, 2.98151566067, -6.31914477675, 9.93943649152, -3.85656214776, 
                0.335403140994, -16.4081976158, 44.9710976818, -62.4405204452, 1.68615917053, 
                0.240182508566, 60.7421927925, -161.122225846, 227.111753543, -0.200930270689, 
                0.169666824743, -225.162964522, 601.035648282, -844.866443191, 0.229361912231, 
                0.118550192463, 840.928416263, -2242.06046986, 3153.12012424, 0.0580826217659, 
                0.0821324054059, -225.211833175, 600.849084896, -844.794584513, 0.0670276007055, 
                0.0565001859137, 60.4236194058, -160.936520965, 226.418721055, 0.039446975412, 
                0.0383668509393, -16.1373305126, 43.1631351373, -60.6285330032, 0.0289844976466, 
                0.026012410329, 4.35985761233, -11.5383128312, 16.2736904204, 0.0198150340018, 
                0.0173635077448, -1.14411787221, 3.10981068213, -4.3400309364, 0.0141553663462, 
                0.0115735586918, 0.322674650716, -0.820795295186, 1.17708400622, 0.00984350061329, 
                0.00750225748812, -0.0759040854886, 0.22829277818, -0.305643224685, 0.00659063120061, 
                0.00476741135574, 0.027192478335, -0.0557235831243, 0.0884888996886, 0.00438397458426, 
                0.00298409708892, -0.00290776333517, 0.0184507156078, -0.0193655210226, 0.00286147046237, 
                0.00182420028857, 0.00351438145734, -0.00274642984457, 0.00794884748406, 0.00181014356628, 
                0.00110710175681, 0.000736547183215, 0.00209184463068, -0.000413152067803, 0.00114195527251, 
                0.000653392684187, 0.000811236261408, 0.000243567450861, 0.00116386831403, 0.000702035343672, 
                0.000380927506439, 0.000382985190506, 0.000408545780925, 0.000330391356467, 0.000431103352798, 
                0.000166220551263, 0.000207040733211, 8.91590990565e-05, 0.000284325397627, 0.000202633875912, 
                0.000108714724369, 0.000110440117265, 0.000123585549038, 9.56548564435e-05, 0.000130720562044, 
                7.89e-05, 8.82580645161e-05, 6.93451612903e-05, 9.77371863799e-05, 9.545e-05 ],
                3,
                3 ]
        self.tkcZonalMin = [ [ 0.0, 0.0, 0.0, 15000.0, 25000.0, 
                35000.0, 45000.0, 55000.0, 65000.0, 80000.0, 
                80000.0, 80000.0 ],
                [ -1.0, -8.65347787874, 13.1081150504, 21.6524689906, 24.9770710059, 
                36.4851049739, 52.1122991508, 70.5233003639, 29.0, 0.0, 
                0.0, 0.0 ],
                2 ]
        self.tkcZonalMean = [ [ 0.0, 0.0, 0.0, 15000.0, 25000.0, 
                35000.0, 45000.0, 55000.0, 65000.0, 80000.0, 
                80000.0, 80000.0 ],
                [ -2.20686521207e-16, -5.89525513109, 28.7555953059, 34.4315132084, 40.6553254438, 
                57.6365341289, 81.5254697829, 101.131941522, 66.0, 0.0, 
                0.0, 0.0 ],
                2 ]
        self.tkcZonalMax = [ [ 0.0, 0.0, 0.0, 15000.0, 25000.0, 
                35000.0, 45000.0, 55000.0, 65000.0, 80000.0, 
                80000.0, 80000.0 ],
                [ 1.0, -3.13703238345, 44.4030755614, 47.2105574262, 56.3335798817, 
                78.7879632839, 110.938640415, 131.740582679, 103.0, 0.0, 
                0.0, 0.0 ],
                2 ]
        self.tkcMeridionalMin = [ [ 0.0, 0.0, 0.0, 15000.0, 25000.0, 
                35000.0, 45000.0, 55000.0, 65000.0, 80000.0, 
                80000.0, 80000.0 ],
                [ 0.384615384615, 3.32826072259, -5.04158271169, -8.3278726887, -9.60656577151, 
                -14.0327326823, -20.0431919811, -27.1243462938, -11.1538461538, 0.0, 
                0.0, 0.0 ],
                2 ]
        self.tkcMeridionalMean = [ [ 0.0, 0.0, 0.0, 15000.0, 25000.0, 
                35000.0, 45000.0, 55000.0, 65000.0, 80000.0, 
                80000.0, 80000.0 ],
                [ 0.284050811574, 1.3763292252, 1.55288530117, 0.167307456216, 0.364935452013, 
                0.211887171154, 0.156260928264, -1.47588115798, 3.87261820748, 0.0, 
                0.0, 0.0 ],
                2 ]
        self.tkcMeridionalMax = [ [ 0.0, 0.0, 0.0, 15000.0, 25000.0, 
                35000.0, 45000.0, 55000.0, 65000.0, 80000.0, 
                80000.0, 80000.0 ],
                [ 0.183486238532, -0.575602272193, 8.14735331402, 8.66248760113, 10.3364366755, 
                14.4565070246, 20.3557138376, 24.1725839779, 18.8990825688, 0.0, 
                0.0, 0.0 ],
                2 ]
        self.tkcUncertaintyTemp = [ [ 0.0, 0.0, 10000.0, 55000.0, 75000.0, 
                100000.0, 100000.0 ],
                [ 5.0, 10.0, 20.0, 20.0, 10.0, 
                0.0, 0.0 ],
                1 ]
        self.tkcUncertaintyPres = [ [ 0.0, 0.0, 33000.0, 60000.0, 90000.0, 
                150000.0, 150000.0 ],
                [ 0.02, 0.1, 0.3, 0.3, 0.45, 
                0.0, 0.0 ],
                1 ]
        self.tkcUncertaintyDens = [ [ 0.0, 0.0, 33000.0, 60000.0, 90000.0, 
                150000.0, 150000.0 ],
                [ 0.02, 0.1, 0.3, 0.3, 0.45, 
                0.0, 0.0 ],
                1 ]