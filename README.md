# ParametricAcousticArray
MATLAB routines for solving the Parametric Acoustic Array


## What's included
Within the download you'll find the following directories and files. 
You'll see something like this:
```
├── README.md
├── absorp
│   ├── cal_absorp_coeff.m
│   └── fig
│       ├── absorp_coeff.jpg
│       └── absorp_coeff.m
├── direct
│   ├── cal_beamwidth.m
│   ├── cal_direct.m
│   ├── cal_effect_len.m
│   └── fig
│       ├── cmp_beamwidth.m
│       ├── cmp_beamwidth_hr50_Tc25.jpg
│       ├── cmp_direct.m
│       ├── cmp_direct_fA1k_fU60k_a0p1_hr50_Tc25.jpg
│       ├── cmp_direct_fA2k_fU60k_a0p1_hr50_Tc25.jpg
│       ├── cmp_direct_fA4k_fU60k_a0p1_hr50_Tc25.jpg
│       ├── cmp_direct_fA500_fU60k_a0p1_hr50_Tc25.jpg
│       ├── cmp_effect_len.m
│       └── cmp_effect_len_Tc20.jpg
├── kzk
│   ├── data
│   │   ├── kzk_cal_gen_fb1k_f160k_a0p1_P0118.mat
│   │   ├── kzk_cal_gen_fb1k_f160k_a0p1_P0135.mat
│   │   ├── kzk_cal_gen_fb2k_f160k_a0p1_P0118.mat
│   │   ├── kzk_cal_gen_fb2k_f160k_a0p1_P0135.mat
│   │   ├── kzk_cal_gen_fb4k_f160k_a0p1_P0118.mat
│   │   ├── kzk_cal_gen_fb4k_f160k_a0p1_P0135.mat
│   │   ├── kzk_cal_gen_fb500_f160k_a0p1_P0118.mat
│   │   └── kzk_cal_gen_fb500_f160k_a0p1_P0135.mat
│   ├── fig
│   │   ├── kzk_plt1.m
│   │   ├── kzk_plt1_set1_P0118.jpg
│   │   ├── kzk_plt1_set2_P0135.jpg
│   │   ├── kzk_plt2.m
│   │   ├── kzk_plt2_cache.jpg
│   │   ├── kzk_plt2_set1_P0118.jpg
│   │   └── kzk_plt2_set2_P0135.jpg
│   ├── kzk_cal.m
│   ├── kzk_cal_gen.m
│   └── kzk_import_z_to_mat.m
└── sandbox.m
```

## Hints for Use
- Run `sandbox.m` first for initialization.
- Run the m files in the subfolder `/fig/` to generate the specific 
figures.
	- Eg.: Run `kzk/fig/kzk_plt1.m` to generate the following figure:
	![Fig. 1](https://github.com/JiaxinZhong/ParametricAcousticArray/blob/master/kzk/fig/kzk_plt1_set2_P0135.jpg)
