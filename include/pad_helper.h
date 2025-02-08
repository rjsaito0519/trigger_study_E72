#ifndef PAD_HELPER_HH
#define PAD_HELPER_HH

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <TMath.h>
#include <TVector3.h>

namespace pad_helper
{
    static const std::vector<Int_t> noisy_pad{ // 1 origin
        0, 50, 56, 57, 58, 59, 68, 72, 73, 74, 76, 110, 108, 109, 137, 139, 183, 224, 226, 282,
        337, 335, 472, 470, 552, 553, 631, 629, 723, 651, 814, 812, 854, 855, 856, 857, 858,
        859, 860, 861, 862, 863, 864, 918, 919, 1042, 1044, 1058, 1059, 1060, 1061, 1062, 1063,
        1094, 1137, 1138, 1288, 1289, 1290, 1291, 1361, 1362, 1363, 1444, 1446, 1448, 1455, 1509,
        1510, 1511, 1512, 1513, 1519, 1523, 1531, 1532, 1543, 1558, 1562, 1563, 1564, 1565, 1584, 
        1586, 1587, 1588, 1589, 1590, 1591, 1592, 1593, 1594, 1595, 1608, 1610, 1612, 1619, 1620, 
        1626, 1630, 1632, 1634, 1636, 1638, 1640, 1642, 1644, 1646, 1662, 1678, 1701, 1705, 1706, 
        1714, 1715, 1733, 1734, 1735, 1743, 1751, 1755, 1757, 1767, 1777, 1778, 1820, 1852, 1853, 
        2011, 2034, 2035, 2036, 2037, 2038, 2039, 2051, 2055, 2059, 2071, 2075, 2217, 2219, 2242, 
        2243, 2244, 2245, 2246, 2247, 2248, 2278, 2280, 2297, 2299, 2303, 2453, 2454, 2455, 2456, 
        2457, 2467, 2469, 2470, 2471, 2472, 2483, 2485, 2487, 2489, 2667, 2669, 2670, 2671, 2675, 
        2681, 2701, 2702, 2703, 2790, 2791, 2792, 2793, 2794, 2795, 2796, 2797, 2798, 2799, 2800, 
        2801, 2802, 2803, 2804, 2805, 2806, 2807, 2808, 2809, 2810, 2811, 2812, 2813, 2814, 2815, 
        2816, 2846, 2885, 2887, 2888, 2889, 2919, 2921, 2999, 3000, 3001, 3002, 3003, 3004, 3005, 
        3006, 3007, 3008, 3009, 3010, 3011, 3012, 3013, 3014, 3015, 3045, 3046, 3047, 3048, 3049, 
        3050, 3051, 3052, 3053, 3054, 3055, 3056, 3057, 3058, 3110, 3111, 3112, 3113, 3143, 3145, 
        3227, 3228, 3229, 3230, 3291, 3292, 3293, 3294, 3295, 3296, 3297, 3298, 3299, 3300, 3301, 
        3314, 3315, 3341, 3342, 3343, 3344, 3345, 3346, 3374, 3375, 3376, 3377, 3414, 3415, 3416, 
        3417, 3418, 3419, 3420, 3421, 3422, 3449, 3451, 3452, 3467, 3473, 3475, 3478, 3479, 3480, 
        3481, 3482, 3483, 3525, 3526, 3527, 3528, 3529, 3530, 3541, 3542, 3543, 3544, 3545, 3546, 
        3547, 3548, 3549, 3565, 3566, 3567, 3568, 3569, 3570, 3572, 3578, 3580, 3581, 3605, 3609, 
        3612, 3614, 3661, 3662, 3663, 3664, 3665, 3666, 3667, 3668, 3669, 3720, 3722, 3724, 3737, 
        3739, 3783, 3784, 3785, 3787, 3788, 3796, 3797, 3798, 3799, 3814, 3815, 3816, 3846, 3848, 
        3905, 3906, 3907, 3908, 3909, 3910, 3911, 4034, 4035, 4036, 4037, 4067, 4069, 4105, 4106, 
        4107, 4108, 4135, 4136, 4137, 4138, 4139, 4140, 4141, 4171, 4172, 4173, 4174, 4277, 4278, 
        4279, 4280, 4354, 4355, 4356, 4357, 4358, 4359, 4385, 4386, 4387, 4415, 4416, 4417, 4418, 
        4419, 4420, 4421, 4422, 4423, 4424, 4425, 4426, 4427, 4428, 4429, 4430, 4431, 4432, 4433, 
        4434, 4435, 4436, 4437, 4438, 4439, 4440, 4441, 4442, 4443, 4444, 4445, 4446, 4447, 4448, 
        4449, 4482, 4484, 4492, 4496, 4569, 4570, 4571, 4572, 4600, 4611, 4612, 4613, 4614, 4615, 
        4616, 4617, 4618, 4619, 4620, 4621, 4622, 4623, 4624, 4625, 4626, 4627, 4663, 4664, 4665, 
        4682, 4684, 4720, 4721, 4722, 4723, 4724, 4725, 4726, 4727, 4728, 4775, 4776, 4777, 4778, 
        4779, 4780, 4879, 4881, 4927, 4928, 4929, 4930, 4931, 4932, 4933, 4934, 4935, 4936, 4937, 
        4945, 4947, 4981, 4982, 4983, 4984, 4985, 5035, 5073, 5074, 5075, 5076, 5092, 5098, 5102, 
        5103, 5104, 5112, 5113, 5114, 5115, 5132, 5133, 5134, 5135, 5136, 5137, 5138, 5139, 5140, 
        5141, 5183, 5184, 5185, 5186, 5187, 5214, 5215, 5216, 5217, 5297, 5298, 5299, 5328, 5329, 
        5330, 5332, 5375, 5376, 5377, 5378, 5406, 5407, 5445, 5447, 5448, 5449, 5450, 5451, 5488, 
        5489, 5490, 5491, 5492, 5533, 5534, 5535, 5562, 5563, 5564, 5565, 5566, 5611, 5612, 5613, 
        5614, 5615, 5616, 5655, 5656, 5657, 5658, 5659, 5668, 5715, 5716, 5717, 5718, 5719, 5720, 
        5721, 5740, 5741, 5758, 5759, 5760, 5761, 5762, 5766, 5767, 5768
    };

    static const Int_t NumOfLayersTPC = 32;
    static const Int_t NoiseChannel[NumOfLayersTPC][42] =
    {
        {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},//0
        {8,9,10,20,24,25,26,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {12,13,41,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {15,56,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {18,71,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {21,86,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //5
        {24,25,101,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {27,116,158,159,160,161,162,163,164,165,166,167,168,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {30,31,154,170,171,172,173,174,175,206,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {33,34,184,185,186,187,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {17,40,42,43,44,46,165,166,167,168,169,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //10
        {10,11,34,35,36,37,38,40,56,78,154,177,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {6,80,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {34,35,37,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {3,28,29,30,31,32,64,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {27,28,29,41,43,44,57,61,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //15
        {27,61,150,152,154,155,156,157,158,159,160,161,162,164,165,166,167,168,169,170,171,172,173,174,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {25,27,59,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,185,190,191,192,193,194,195,196,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {2,26,27,59,143,144,207,208,209,210,211,212,213,214,215,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {25,26,28,58,59,98,99,100,101,102,103,104,133,135,136,151,157,159,162,163,164,165,166,209,210,211,212,225,226,227,228,229,230,231,-1,-1,-1,-1,-1,-1,-1,-1},
        {11,12,14,15,16,24,26,58,107,108,109,110,111,112,113,114,242,243,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //20
        {16,48,107,108,109,110,111,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {4,5,37,75,76,105,106,107,108,109,141,142,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {29,30,106,107,108,109,137,138,139,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199},
        {24,153,154,155,156,157,158,159,160,161,162,163,164,165,166,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {18,56,57,58,59,60,61,62,111,112,113,114,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //25
        {13,61,62,63,64,65,66,68,69,79,115,116,117,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {7,8,36,46,47,66,67,68,69,70,71,72,73,117,118,119,148,149,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {35,113,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {5,8,9,48,93,122,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {2,41,42,43,44,85,86,87,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //30
        {37,38,39,40,62,63,80,81,82,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}
    };

    //#PadID is defined as 0 origin
    //#OfPad #division #radius padLength
    static const Double_t padParameter[32][6] =
    {
        {0, 48,    14.75, 48, 0,  9.},
        {1, 48,    24.25, 48, 0,  9.},
        {2, 72,    33.75, 72, 0,  9.},
        {3, 96,    43.25, 96, 0,  9.},
        {4, 120,    52.75,120,0,   9.},
        {5, 144,    62.25,144,0,   9.},
        {6, 168,    71.75,168,0,   9.},
        {7, 192,    81.25,192,0,   9.},
        {8, 216,    90.75,216,0,   9.},
        {9, 240,    100.25,240,0,  9.},
        {10,208,    111.5,241, 0,  12.5},
        {11,218,    124.5,271, 0,  12.5},
        {12,230,    137.5,300, 0,  12.5},
        {13,214,    150.5,330, 0,  12.5},
        {14,212,    163.5,360, 0,  12.5},
        {15,214,    176.5,390, 0,  12.5},
        {16,220,    189.5,420, 0,  12.5},
        {17,224,    202.5,449, 0,  12.5},
        {18,232,    215.5,479, 0,  12.5},
        {19,238,    228.5,509, 0,  12.5},
        {20,244,    241.5,539, 0,  12.5},
        {21,232,    254.5,569, 0,  12.5},
        {22,218,    267.5,599, 0,  12.5},
        {23,210,    280.5,628, 0,  12.5},
        {24,206,    293.5,658, 0,  12.5},
        {25,202,    306.5,688, 0,  12.5},
        {26,200,    319.5,718, 0,  12.5},
        {27,196,    332.5,748, 0,  12.5},
        {28,178,    345.5,777, 0,  12.5},
        {29,130,    358.5,807, 0,  12.5},
        {30,108,    371.5,837, 0,  12.5},
        {31,90,     384.5,867, 0, 12.5}
    };

    inline Double_t getDTheta(Int_t layerID)
    {
        return (360. / padParameter[layerID][3]);
    }

    inline Double_t getsTheta(Int_t layerID)
    {
        Double_t sTheta = 180. - (360. / padParameter[layerID][3]) * padParameter[layerID][1] / 2.;
        return sTheta;
    }

    inline Double_t getRadius(Int_t layerID)
    {
        return padParameter[layerID][2];
    }

    inline Double_t getLength(Int_t layerID)
    {
        return padParameter[layerID][5];
    }

    inline Int_t getPadID(Int_t layerID, Int_t rowID)
    {
        Int_t padID = 0;
        for (int layi = 0; layi < layerID; layi++) padID += padParameter[layi][1];
        padID += rowID;
        return padID;
    }

    inline Int_t getLayerID(Int_t padID)
    {
        //    padID-=1;
        int layer;
        int sum = 0;

        for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
        {
            sum += padParameter[layer][1];
        }
        return layer;
    }

    inline Int_t getRowID(Int_t padID)
    {
        //    padID-=1;
        int layer, row;
        int sum = 0;

        for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
        {
            sum += padParameter[layer][1];
        }
        row = padID - sum;
        return row;
    }

    inline int findPadID(double z, double x)
    {
        z += 143;
        double radius = sqrt(x * x + z * z);
        double angle;
        if (z == 0)
        {
            if (x > 0)   angle = 1.5 * TMath::Pi();
            else if (x < 0)   angle = 0.5 * TMath::Pi();
            else return -1000; // no padID if (0,0)
        }
        else
        {
            if (z > 0) angle = TMath::Pi() + atan(x / z);
            else if (z < 0 && x < 0) angle = atan(x / z);
            else angle = 2 * TMath::Pi() + atan(x / z);
            //	angle = TMath::Pi() - atan(-x / z);
        }

        int layer, row;
        // find layer_num.
        for (layer = 0; !(padParameter[layer][2] + padParameter[layer][5] * 0.5 >= radius
            && padParameter[layer][2] - padParameter[layer][5] * 0.5 <= radius); layer++)
        {
            if (layer >= 32) return -1000;
            if (layer != 0)
            {
                if (padParameter[layer][2] - padParameter[layer][5] * 0.5 >= radius &&
                    padParameter[layer - 1][2] + padParameter[layer - 1][5] * 0.5 <= radius) return -layer;
            }
        }

        //std::cout<<"padHelper:: layer="<<layer<<", angle="<<angle<<", "<<(getsTheta(layer)*TMath::Pi()/180.)<<std::endl;
        // find row_num
        //  if (angle - (padParameter[layer][4]*TMath::Pi()/180.) < 0) return -1000;
        if (angle - (getsTheta(layer) * TMath::Pi() / 180.) < 0) return -2000;

        //    double a, b, c;
        //row = (int)((angle-(padParameter[layer][4]*TMath::Pi()/180.))/(padParameter[layer][3]*TMath::Pi()/180.));
        //    row = (int)((angle-(getsTheta(layer)*TMath::Pi()/180.))/(getDTheta(layer)*TMath::Pi()/180.));

        //row = (int)((angle-(getsTheta(layer)*TMath::Pi()/180.))/(getDTheta(layer)*TMath::Pi()/180.))+1;
        row = (int)((angle - (getsTheta(layer) * TMath::Pi() / 180.)) / (getDTheta(layer) * TMath::Pi() / 180.));
        if (row > padParameter[layer][1]) return -1000;

        return getPadID(layer, row);
    }

    inline Double_t getTheta(Int_t padID)
    {
        //    padID-=1;
        int layer, row;
        int sum = 0;

        for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
        {
            sum += padParameter[layer][1];
        }
        row = padID - sum;
        Double_t sTheta = 180. - (360. / padParameter[layer][3]) * padParameter[layer][1] / 2.;
        Double_t theta = sTheta + (row + 0.5) * (360. - 2 * sTheta) / padParameter[layer][1];
        return theta;
    }

    inline Double_t getR(Int_t padID)
    {
        //    padID-=1;
        int layer;//, row;
        int sum = 0;

        for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
        {
            sum += padParameter[layer][1];
        }
        //row = padID - sum;
        Double_t R = padParameter[layer][2];
        return R;
    }

    inline TVector3 getPoint(int padID)
    {
        // 0 originでOK, 1 originなら-1を入れる
        // padID-=1;
        int layer, row;
        int sum = 0;

        for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
        {
            sum += padParameter[layer][1];
        }
        row = padID - sum;

        TVector3 result;
        if (row > padParameter[layer][1])
        { // out of range
            result.SetX(0);
            result.SetY(-1);
            result.SetZ(0);
        }
        else
        {
            double x, z;
            Double_t sTheta = 180. - (360. / padParameter[layer][3]) * padParameter[layer][1] / 2.;
            x = padParameter[layer][2] * -sin((360. / padParameter[layer][3]) * TMath::Pi() / 180. * (row + 0.5) + sTheta * TMath::Pi() / 180.);
            z = padParameter[layer][2] * -cos((360. / padParameter[layer][3]) * TMath::Pi() / 180. * (row + 0.5) + sTheta * TMath::Pi() / 180.) - 143.0;
            result.SetX(x);
            result.SetY(0);
            result.SetZ(z);
        }
        return result;
    }
}

#endif