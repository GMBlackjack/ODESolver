# As documented in the NRPy+ tutorial module
#   Tutorial-RK_Butcher_Table_Dictionary.ipynb ,
#   this module will construct a Python dictionary
#   of Butcher tables for a large number of explicit
#   RK methods.

# Authors: Brandon Clark
#          Zachariah B. Etienne
#          zachetie **at** gmail **dot* com

# Step 1: Initialize needed Python/NRPy+ modules
import sympy as sp  # Import SymPy, a computer algebra system written entirely in Python

# Step 2a: Generating a Dictionary of Butcher Tables for Explicit Runge Kutta Techniques

# Initialize the dictionary Butcher_dict
Butcher_dict = {}

# Step 2.a.i: Euler's Method

Butcher_dict['Euler'] = (
[[sp.sympify(0)],
 ["", sp.sympify(1)]]
, 1)

# Step 2.a.ii: RK2 Heun's Method

Butcher_dict['RK2 Heun'] = (
[[sp.sympify(0)],
[sp.sympify(1), sp.sympify(1)],
["", sp.Rational(1,2), sp.Rational(1,2)]]
, 2)

# Step 2.a.iii: RK2 Midpoint (MP) Method

Butcher_dict['RK2 MP'] = (
[[sp.sympify(0)],
[sp.Rational(1,2), sp.Rational(1,2)],
["", sp.sympify(0), sp.sympify(1)]]
, 2)

# Step 2.a.iv: RK2 Ralston's Method

Butcher_dict['RK2 Ralston'] = (
[[sp.sympify(0)],
[sp.Rational(2,3), sp.Rational(2,3)],
["", sp.Rational(1,4), sp.Rational(3,4)]]
, 2)

# Step 2.a.v: Kutta's  Third-order Method

Butcher_dict['RK3'] = (
[[sp.sympify(0)],
[sp.Rational(1,2), sp.Rational(1,2)],
[sp.sympify(1), sp.sympify(-1), sp.sympify(2)],
["", sp.Rational(1,6), sp.Rational(2,3), sp.Rational(1,6)]]
, 3)

# Step 2.a.vi: RK3 Heun's Method

Butcher_dict['RK3 Heun'] = (
[[sp.sympify(0)],
[sp.Rational(1,3), sp.Rational(1,3)],
[sp.Rational(2,3), sp.sympify(0), sp.Rational(2,3)],
["", sp.Rational(1,4), sp.sympify(0), sp.Rational(3,4)]]
, 3)

# Step 2.a.vii: RK3 Ralton's Method

Butcher_dict['RK3 Ralston'] = ([[0],
[sp.Rational(1,2), sp.Rational(1,2)],
[sp.Rational(3,4), sp.sympify(0), sp.Rational(3,4)],
["", sp.Rational(2,9), sp.Rational(1,3), sp.Rational(4,9)]]
, 3)

# Step 2.a.viii: Strong Stability Preserving Runge-Kutta (SSPRK3) Method
Butcher_dict['SSPRK3'] = (
[[0],
[sp.sympify(1), sp.sympify(1)],
[sp.Rational(1,2), sp.Rational(1,4), sp.Rational(1,4)],
["", sp.Rational(1,6), sp.Rational(1,6), sp.Rational(2,3)]]
, 3)

# Step 2.a.ix: Classic RK4 Method

Butcher_dict['RK4'] = ([[sp.sympify(0)],
[sp.Rational(1,2), sp.Rational(1,2)],
[sp.Rational(1,2), sp.sympify(0), sp.Rational(1,2)],
[sp.sympify(1), sp.sympify(0), sp.sympify(0), sp.sympify(1)],
["", sp.Rational(1,6), sp.Rational(1,3), sp.Rational(1,3), sp.Rational(1,6)]]
, 4)

# Step 2.a.x:  RK5 Dormand-Prince Method

Butcher_dict['DP5'] = ([[0],
[sp.Rational(1,5), sp.Rational(1,5)],
[sp.Rational(3,10),sp.Rational(3,40), sp.Rational(9,40)],
[sp.Rational(4,5), sp.Rational(44,45), sp.Rational(-56,15), sp.Rational(32,9)],
[sp.Rational(8,9), sp.Rational(19372,6561), sp.Rational(-25360,2187), sp.Rational(64448,6561), sp.Rational(-212,729)],
[sp.sympify(1), sp.Rational(9017,3168), sp.Rational(-355,33), sp.Rational(46732,5247), sp.Rational(49,176), sp.Rational(-5103,18656)],
[sp.sympify(1), sp.Rational(35,384), sp.sympify(0), sp.Rational(500,1113), sp.Rational(125,192), sp.Rational(-2187,6784), sp.Rational(11,84)],
["", sp.Rational(35,384), sp.sympify(0), sp.Rational(500,1113), sp.Rational(125,192), sp.Rational(-2187,6784), sp.Rational(11,84), sp.sympify(0)]]
, 5)

# Step 2.a.xi:  RK5 Dormand-Prince Method Alternative

Butcher_dict['DP5alt'] = (
[[0],
[sp.Rational(1,10), sp.Rational(1,10)],
[sp.Rational(2,9), sp.Rational(-2, 81), sp.Rational(20, 81)],
[sp.Rational(3,7), sp.Rational(615, 1372), sp.Rational(-270, 343), sp.Rational(1053, 1372)],
[sp.Rational(3,5), sp.Rational(3243, 5500), sp.Rational(-54, 55), sp.Rational(50949, 71500), sp.Rational(4998, 17875)],
[sp.Rational(4, 5), sp.Rational(-26492, 37125), sp.Rational(72, 55), sp.Rational(2808, 23375), sp.Rational(-24206, 37125), sp.Rational(338, 459)],
[sp.sympify(1), sp.Rational(5561, 2376), sp.Rational(-35, 11), sp.Rational(-24117, 31603), sp.Rational(899983, 200772), sp.Rational(-5225, 1836), sp.Rational(3925, 4056)],
["", sp.Rational(821, 10800), sp.sympify(0), sp.Rational(19683, 71825), sp.Rational(175273, 912600), sp.Rational(395, 3672), sp.Rational(785, 2704), sp.Rational(3, 50)]]
, 5)

# Step 2.a.xii:  RK5 Dormand-Prince Method Alternative
Butcher_dict['DP5alt'] = (
[[0],
[sp.Rational(1,10), sp.Rational(1,10)],
[sp.Rational(2,9), sp.Rational(-2, 81), sp.Rational(20, 81)],
[sp.Rational(3,7), sp.Rational(615, 1372), sp.Rational(-270, 343), sp.Rational(1053, 1372)],
[sp.Rational(3,5), sp.Rational(3243, 5500), sp.Rational(-54, 55), sp.Rational(50949, 71500), sp.Rational(4998, 17875)],
[sp.Rational(4, 5), sp.Rational(-26492, 37125), sp.Rational(72, 55), sp.Rational(2808, 23375), sp.Rational(-24206, 37125), sp.Rational(338, 459)],
[sp.sympify(1), sp.Rational(5561, 2376), sp.Rational(-35, 11), sp.Rational(-24117, 31603), sp.Rational(899983, 200772), sp.Rational(-5225, 1836), sp.Rational(3925, 4056)],
["", sp.Rational(821, 10800), sp.sympify(0), sp.Rational(19683, 71825), sp.Rational(175273, 912600), sp.Rational(395, 3672), sp.Rational(785, 2704), sp.Rational(3, 50)]]
, 5)

# Step 2.a.xiii:  RK5 Cash-Karp Method

Butcher_dict['CK5'] = (
[[0],
[sp.Rational(1,5), sp.Rational(1,5)],
[sp.Rational(3,10),sp.Rational(3,40), sp.Rational(9,40)],
[sp.Rational(3,5), sp.Rational(3,10), sp.Rational(-9,10), sp.Rational(6,5)],
[sp.sympify(1), sp.Rational(-11,54), sp.Rational(5,2), sp.Rational(-70,27), sp.Rational(35,27)],
[sp.Rational(7,8), sp.Rational(1631,55296), sp.Rational(175,512), sp.Rational(575,13824), sp.Rational(44275,110592), sp.Rational(253,4096)],
["",sp.Rational(37,378), sp.sympify(0), sp.Rational(250,621), sp.Rational(125,594), sp.sympify(0), sp.Rational(512,1771)]]
, 5)

# Step 2.a.xiv:  RK6 Dormand-Prince Method

Butcher_dict['DP6'] = (
[[0],
[sp.Rational(1,10), sp.Rational(1,10)],
[sp.Rational(2,9), sp.Rational(-2, 81), sp.Rational(20, 81)],
[sp.Rational(3,7), sp.Rational(615, 1372), sp.Rational(-270, 343), sp.Rational(1053, 1372)],
[sp.Rational(3,5), sp.Rational(3243, 5500), sp.Rational(-54, 55), sp.Rational(50949, 71500), sp.Rational(4998, 17875)],
[sp.Rational(4, 5), sp.Rational(-26492, 37125), sp.Rational(72, 55), sp.Rational(2808, 23375), sp.Rational(-24206, 37125), sp.Rational(338, 459)],
[sp.sympify(1), sp.Rational(5561, 2376), sp.Rational(-35, 11), sp.Rational(-24117, 31603), sp.Rational(899983, 200772), sp.Rational(-5225, 1836), sp.Rational(3925, 4056)],
[sp.sympify(1), sp.Rational(465467, 266112), sp.Rational(-2945, 1232), sp.Rational(-5610201, 14158144), sp.Rational(10513573, 3212352), sp.Rational(-424325, 205632), sp.Rational(376225, 454272), sp.sympify(0)],
["", sp.Rational(61, 864), sp.sympify(0), sp.Rational(98415, 321776), sp.Rational(16807, 146016), sp.Rational(1375, 7344), sp.Rational(1375, 5408), sp.Rational(-37, 1120), sp.Rational(1,10)]]
, 6)

#  Step 2.a.xv:  RK6 Luther's Method

q = sp.sqrt(21)
Butcher_dict['L6'] = (
[[0],
[sp.sympify(1), sp.sympify(1)],
[sp.Rational(1,2), sp.Rational(3,8), sp.Rational(1,8)],
[sp.Rational(2,3), sp.Rational(8,27), sp.Rational(2,27), sp.Rational(8,27)],
[(7 - q)/14, (-21 + 9*q)/392, (-56 + 8*q)/392, (336 -48*q)/392, (-63 + 3*q)/392],
[(7 + q)/14, (-1155 - 255*q)/1960, (-280 -  40*q)/1960, (-320*q)/1960, (63 + 363*q)/1960, (2352 + 392*q)/1960],
[sp.sympify(1), ( 330 + 105*q)/180, sp.Rational(2,3), (-200 + 280*q)/180, (126 - 189*q)/180, (-686 - 126*q)/180, (490 -  70*q)/180],
["", sp.Rational(1, 20), sp.sympify(0), sp.Rational(16, 45), sp.sympify(0), sp.Rational(49, 180), sp.Rational(49, 180), sp.Rational(1, 20)]]
, 6)

# Step 2.a.xvi: RK8 Dormand-Prince Method

Butcher_dict['DP8']=(
[[0],
[sp.Rational(1, 18), sp.Rational(1, 18)],
[sp.Rational(1, 12), sp.Rational(1, 48), sp.Rational(1, 16)],
[sp.Rational(1, 8), sp.Rational(1, 32), sp.sympify(0), sp.Rational(3, 32)],
[sp.Rational(5, 16), sp.Rational(5, 16), sp.sympify(0), sp.Rational(-75, 64), sp.Rational(75, 64)],
[sp.Rational(3, 8), sp.Rational(3, 80), sp.sympify(0), sp.sympify(0), sp.Rational(3, 16), sp.Rational(3, 20)],
[sp.Rational(59, 400), sp.Rational(29443841, 614563906), sp.sympify(0), sp.sympify(0), sp.Rational(77736538, 692538347), sp.Rational(-28693883, 1125000000), sp.Rational(23124283, 1800000000)],
[sp.Rational(93, 200), sp.Rational(16016141, 946692911), sp.sympify(0), sp.sympify(0), sp.Rational(61564180, 158732637), sp.Rational(22789713, 633445777), sp.Rational(545815736, 2771057229), sp.Rational(-180193667, 1043307555)],
[sp.Rational(5490023248, 9719169821), sp.Rational(39632708, 573591083), sp.sympify(0), sp.sympify(0), sp.Rational(-433636366, 683701615), sp.Rational(-421739975, 2616292301), sp.Rational(100302831, 723423059), sp.Rational(790204164, 839813087), sp.Rational(800635310, 3783071287)],
[sp.Rational(13, 20), sp.Rational(246121993, 1340847787), sp.sympify(0), sp.sympify(0), sp.Rational(-37695042795, 15268766246), sp.Rational(-309121744, 1061227803), sp.Rational(-12992083, 490766935), sp.Rational(6005943493, 2108947869), sp.Rational(393006217, 1396673457), sp.Rational(123872331, 1001029789)],
[sp.Rational(1201146811, 1299019798), sp.Rational(-1028468189, 846180014), sp.sympify(0), sp.sympify(0), sp.Rational(8478235783, 508512852), sp.Rational(1311729495, 1432422823), sp.Rational(-10304129995, 1701304382), sp.Rational(-48777925059, 3047939560), sp.Rational(15336726248, 1032824649), sp.Rational(-45442868181, 3398467696), sp.Rational(3065993473, 597172653)],
[sp.sympify(1), sp.Rational(185892177, 718116043), sp.sympify(0), sp.sympify(0), sp.Rational(-3185094517, 667107341), sp.Rational(-477755414, 1098053517), sp.Rational(-703635378, 230739211), sp.Rational(5731566787, 1027545527), sp.Rational(5232866602, 850066563), sp.Rational(-4093664535, 808688257), sp.Rational(3962137247, 1805957418), sp.Rational(65686358, 487910083)],
[sp.sympify(1), sp.Rational(403863854, 491063109), sp.sympify(0), sp.sympify(0), sp.Rational(-5068492393, 434740067), sp.Rational(-411421997, 543043805), sp.Rational(652783627, 914296604), sp.Rational(11173962825, 925320556), sp.Rational(-13158990841, 6184727034), sp.Rational(3936647629, 1978049680), sp.Rational(-160528059, 685178525), sp.Rational(248638103, 1413531060), sp.sympify(0)],
["", sp.Rational(14005451, 335480064), sp.sympify(0), sp.sympify(0), sp.sympify(0), sp.sympify(0), sp.Rational(-59238493, 1068277825), sp.Rational(181606767, 758867731), sp.Rational(561292985, 797845732), sp.Rational(-1041891430, 1371343529), sp.Rational(760417239, 1151165299), sp.Rational(118820643, 751138087), sp.Rational(-528747749, 2220607170), sp.Rational(1, 4)]]
, 8)

# TESTING for ADAPTIVE METHODS

Butcher_dict['AHE'] = (
[[sp.sympify(0)],
[sp.sympify(1), sp.sympify(1)],
["", sp.Rational(1,2), sp.Rational(1,2)],
["", sp.sympify(1), sp.sympify(0)]]
, 2)

Butcher_dict['ABS'] = (
[[sp.sympify(0)],
[sp.Rational(1,2), sp.Rational(1,2)],
[sp.Rational(3,4), sp.sympify(0), sp.Rational(3,4)],
[sp.sympify(1), sp.Rational(2,9), sp.Rational(1,3), sp.Rational(4,9)],
["", sp.Rational(2,9), sp.Rational(1,3), sp.Rational(4,9), sp.sympify(0)],
["", sp.Rational(7,24), sp.Rational(1,4), sp.Rational(1,3), sp.Rational(1,8)]]
, 3)

Butcher_dict['ARKF'] = (
[[sp.sympify(0)],
[sp.Rational(1,4), sp.Rational(1,4)],
[sp.Rational(3,8), sp.Rational(3,32), sp.Rational(9,32)],
[sp.Rational(12,13), sp.Rational(1932,2197), sp.Rational(-7200,2197), sp.Rational(7296,2197)],
[sp.sympify(1), sp.Rational(439,216), sp.sympify(-8), sp.Rational(3680,513), sp.Rational(-845,4104)],
[sp.Rational(1,2), sp.Rational(-8,27), sp.sympify(2), sp.Rational(-3544,2565), sp.Rational(1859,4104), sp.Rational(-11,40)],
["", sp.Rational(16,135), sp.sympify(0), sp.Rational(6656,12825), sp.Rational(28561,56430), sp.Rational(-9,50), sp.Rational(2,55)],
["", sp.Rational(25,216), sp.sympify(0), sp.Rational(1408,2565), sp.Rational(2197,4104), sp.Rational(-1,5), sp.sympify(0)]]
, 5)

Butcher_dict['ACK'] = (
[[0],
[sp.Rational(1,5), sp.Rational(1,5)],
[sp.Rational(3,10),sp.Rational(3,40), sp.Rational(9,40)],
[sp.Rational(3,5), sp.Rational(3,10), sp.Rational(-9,10), sp.Rational(6,5)],
[sp.sympify(1), sp.Rational(-11,54), sp.Rational(5,2), sp.Rational(-70,27), sp.Rational(35,27)],
[sp.Rational(7,8), sp.Rational(1631,55296), sp.Rational(175,512), sp.Rational(575,13824), sp.Rational(44275,110592), sp.Rational(253,4096)],
["",sp.Rational(37,378), sp.sympify(0), sp.Rational(250,621), sp.Rational(125,594), sp.sympify(0), sp.Rational(512,1771)],
["",sp.Rational(2825,27648), sp.sympify(0), sp.Rational(18575,48384), sp.Rational(13525,55296), sp.Rational(277,14336), sp.Rational(1,4)]]
, 5)

Butcher_dict['ADP5'] = (
[[0],
[sp.Rational(1,5), sp.Rational(1,5)],
[sp.Rational(3,10),sp.Rational(3,40), sp.Rational(9,40)],
[sp.Rational(4,5), sp.Rational(44,45), sp.Rational(-56,15), sp.Rational(32,9)],
[sp.Rational(8,9), sp.Rational(19372,6561), sp.Rational(-25360,2187), sp.Rational(64448,6561), sp.Rational(-212,729)],
[sp.sympify(1), sp.Rational(9017,3168), sp.Rational(-355,33), sp.Rational(46732,5247), sp.Rational(49,176), sp.Rational(-5103,18656)],
[sp.sympify(1), sp.Rational(35,384), sp.sympify(0), sp.Rational(500,1113), sp.Rational(125,192), sp.Rational(-2187,6784), sp.Rational(11,84)],
["", sp.Rational(35,384), sp.sympify(0), sp.Rational(500,1113), sp.Rational(125,192), sp.Rational(-2187,6784), sp.Rational(11,84), sp.sympify(0)],
["", sp.Rational(5179,57600), sp.sympify(0), sp.Rational(7571,16695), sp.Rational(393,640), sp.Rational(-92097,339200), sp.Rational(187,2100), sp.Rational(1,40)]]
, 5)

Butcher_dict['ADP8']=(
[[0],
[sp.Rational(1, 18), sp.Rational(1, 18)],
[sp.Rational(1, 12), sp.Rational(1, 48), sp.Rational(1, 16)],
[sp.Rational(1, 8), sp.Rational(1, 32), sp.sympify(0), sp.Rational(3, 32)],
[sp.Rational(5, 16), sp.Rational(5, 16), sp.sympify(0), sp.Rational(-75, 64), sp.Rational(75, 64)],
[sp.Rational(3, 8), sp.Rational(3, 80), sp.sympify(0), sp.sympify(0), sp.Rational(3, 16), sp.Rational(3, 20)],
[sp.Rational(59, 400), sp.Rational(29443841, 614563906), sp.sympify(0), sp.sympify(0), sp.Rational(77736538, 692538347), sp.Rational(-28693883, 1125000000), sp.Rational(23124283, 1800000000)],
[sp.Rational(93, 200), sp.Rational(16016141, 946692911), sp.sympify(0), sp.sympify(0), sp.Rational(61564180, 158732637), sp.Rational(22789713, 633445777), sp.Rational(545815736, 2771057229), sp.Rational(-180193667, 1043307555)],
[sp.Rational(5490023248, 9719169821), sp.Rational(39632708, 573591083), sp.sympify(0), sp.sympify(0), sp.Rational(-433636366, 683701615), sp.Rational(-421739975, 2616292301), sp.Rational(100302831, 723423059), sp.Rational(790204164, 839813087), sp.Rational(800635310, 3783071287)],
[sp.Rational(13, 20), sp.Rational(246121993, 1340847787), sp.sympify(0), sp.sympify(0), sp.Rational(-37695042795, 15268766246), sp.Rational(-309121744, 1061227803), sp.Rational(-12992083, 490766935), sp.Rational(6005943493, 2108947869), sp.Rational(393006217, 1396673457), sp.Rational(123872331, 1001029789)],
[sp.Rational(1201146811, 1299019798), sp.Rational(-1028468189, 846180014), sp.sympify(0), sp.sympify(0), sp.Rational(8478235783, 508512852), sp.Rational(1311729495, 1432422823), sp.Rational(-10304129995, 1701304382), sp.Rational(-48777925059, 3047939560), sp.Rational(15336726248, 1032824649), sp.Rational(-45442868181, 3398467696), sp.Rational(3065993473, 597172653)],
[sp.sympify(1), sp.Rational(185892177, 718116043), sp.sympify(0), sp.sympify(0), sp.Rational(-3185094517, 667107341), sp.Rational(-477755414, 1098053517), sp.Rational(-703635378, 230739211), sp.Rational(5731566787, 1027545527), sp.Rational(5232866602, 850066563), sp.Rational(-4093664535, 808688257), sp.Rational(3962137247, 1805957418), sp.Rational(65686358, 487910083)],
[sp.sympify(1), sp.Rational(403863854, 491063109), sp.sympify(0), sp.sympify(0), sp.Rational(-5068492393, 434740067), sp.Rational(-411421997, 543043805), sp.Rational(652783627, 914296604), sp.Rational(11173962825, 925320556), sp.Rational(-13158990841, 6184727034), sp.Rational(3936647629, 1978049680), sp.Rational(-160528059, 685178525), sp.Rational(248638103, 1413531060), sp.sympify(0)],
["", sp.Rational(14005451, 335480064), sp.sympify(0), sp.sympify(0), sp.sympify(0), sp.sympify(0), sp.Rational(-59238493, 1068277825), sp.Rational(181606767, 758867731), sp.Rational(561292985, 797845732), sp.Rational(-1041891430, 1371343529), sp.Rational(760417239, 1151165299), sp.Rational(118820643, 751138087), sp.Rational(-528747749, 2220607170), sp.Rational(1, 4)],
["", sp.Rational(13451932, 455176623), sp.sympify(0), sp.sympify(0), sp.sympify(0), sp.sympify(0), sp.Rational(-808719846, 976000145), sp.Rational(1757004468, 5645159321), sp.Rational(656045339, 265891186), sp.Rational(-3867574721, 1518517206), sp.Rational(465885868, 322736535), sp.Rational(53011238, 667516719), sp.Rational(2, 45), sp.sympify(0)]]
, 8)

#Butcher_dict['AB'] = (
#[[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#[3/2, -1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#[23/12, -4/3, 5/12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#[55/24, -59/24, 37/24, -3/8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#[1901/720, -1387/360, 109/30, -637/360, 251/720, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#[4277/1440, -2641/480, 4991/720, -3649/720, 959/480, -95/288, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#[198721/60480, -18637/2520, 235183/20160, -10754/945, 135713/20160, -5603/2520, 19087/60480, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#[16083/4480, -1152169/120960, 242653/13440, -296053/13440, 2102243/120960, -115747/13440, 32863/13440, -5257/17280, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#[14097247/3628800, -21562603/1814400, 47738393/1814400, -69927631/1814400, 862303/22680, -45586321/1814400, 19416743/1814400, -4832053/1814400, 1070017/3628800, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#[4325321/1036800, -104995189/7257600, 6648317/181440, -28416361/453600, 269181919/3628800, -222386081/3628800, 15788639/453600, -2357683/181440, 20884811/7257600, -25713/89600, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#[2132509567/479001600, -2067948781/119750400, 1572737587/31933440, -1921376209/19958400, 3539798831/26611200, -82260679/623700, 2492064913/26611200, -186080291/3991680, 2472634817/159667200, -52841941/17107200, 26842253/95800320, 0, 0, 0, 0, 0, 0, 0, 0], 
#[4527766399/958003200, -6477936721/319334400, 12326645437/191600640, -15064372973/106444800, 35689892561/159667200, -41290273229/159667200, 35183928883/159667200, -625551749/4561920, 923636629/15206400, -17410248271/958003200, 30082309/9123840, -4777223/17418240, 0, 0, 0, 0, 0, 0, 0], 
#[13064406523627/2615348736000, -931781102989/39626496000, 5963794194517/72648576000, -10498491598103/52306974720, 20730767690131/58118860800, -34266367915049/72648576000, 228133014533/486486000, -2826800577631/8072064000, 2253957198793/11623772160, -20232291373837/261534873600, 4588414555201/217945728000, -169639834921/48432384000, 703604254357/2615348736000, 0, 0, 0, 0, 0, 0], 
#[905730205/172204032, -140970750679621/5230697472000, 89541175419277/871782912000, -34412222659093/124540416000, 570885914358161/1046139494400, -31457535950413/38745907200, 134046425652457/145297152000, -350379327127877/435891456000, 310429955875453/581188608000, -10320787460413/38745907200, 7222659159949/74724249600, -21029162113651/871782912000, 6460951197929/1743565824000, -106364763817/402361344000, 0, 0, 0, 0, 0], 
#[13325653738373/2414168064000, -60007679150257/1961511552000, 3966421670215481/31384184832000, -25990262345039/70053984000, 25298910337081429/31384184832000, -2614079370781733/1961511552000, 17823675553313503/10461394944000, -2166615342637/1277025750, 13760072112094753/10461394944000, -1544031478475483/1961511552000, 1600835679073597/4483454976000, -58262613384023/490377888000, 859236476684231/31384184832000, -696561442637/178319232000, 1166309819657/4483454976000, 0, 0, 0, 0], 
#[362555126427073/62768369664000, -2161567671248849/62768369664000, 740161300731949/4828336128000, -4372481980074367/8966909952000, 72558117072259733/62768369664000, -131963191940828581/62768369664000, 62487713370967631/20922789888000, -70006862970773983/20922789888000, 62029181421198881/20922789888000, -129930094104237331/62768369664000, 10103478797549069/8966909952000, -2674355537386529/5706215424000, 9038571752734087/62768369664000, -1934443196892599/62768369664000, 36807182273689/8966909952000, -25221445/98402304, 0, 0, 0], 
#[14845854129333883/2462451425280000, -55994879072429317/1455084933120000, 2612634723678583/14227497123840, -22133884200927593/35177877504000, 5173388005728297701/3201186852864000, -5702855818380878219/1778437140480000, 80207429499737366711/16005934264320000, -3993885936674091251/640237370572800, 2879939505554213/463134672000, -324179886697104913/65330343936000, 7205576917796031023/2286562037760000, -2797406189209536629/1778437140480000, 386778238886497951/640237370572800, -551863998439384493/3201186852864000, 942359269351333/27360571392000, -68846386581756617/16005934264320000, 8092989203533249/32011868528640000, 0, 0], 
#[57424625956493833/9146248151040000, -3947240465864473/92386344960000, 497505713064683651/2286562037760000, -511501877919758129/640237370572800, 65509525475265061/29640619008000, -38023516029116089751/8002967132160000, 129650088885345917773/16005934264320000, -19726972891423175089/1778437140480000, 3146403501110383511/256094948229120, -70617432699294428737/6402373705728000, 14237182892280945743/1778437140480000, -74619315088494380723/16005934264320000, 17195392832483362153/8002967132160000, -4543527303777247/5928123801600, 653581961828485643/3201186852864000, -612172313896136299/16005934264320000, 2460247368070567/547211427840000, -85455477715379/342372925440000, 0], 
#[333374427829017307697/51090942171709440000, -5148905233415267713/109168679854080000, 395276943631267674287/1548210368839680000, -2129159630108649501931/2128789257154560000, 841527158963865085639/283838567620608000, -189774312558599272277/27646613729280000, 856822959645399341657/67580611338240000, -13440468702008745259589/709596419051520000, 196513123964380075325537/8515157028618240000, -57429776853357830333/2494674910728000, 53354279746900330600757/2838385676206080000, -26632588461762447833393/2128789257154560000, 4091553114434184723167/608225502044160000, -291902259907317785203/101370917007360000, 816476630884557765547/851515702861824000, -169944934591213283591/709596419051520000, 239730549209090923561/5676771352412160000, -19963382447193730393/4257578514309120000, 12600467236042756559/51090942171709440000]]
#,19)