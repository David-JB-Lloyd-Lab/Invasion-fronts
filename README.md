# Invasion-fronts
Matlab codes for computing pattern invasion fronts in the planar Swift-Hohenberg equation outside the homoclinic snaking region. 

These codes reproduce the figures in:

Invasion Fronts Outside the Homoclinic Snaking Region in the Planar Swift--Hohenberg Equation
DJB Lloyd, SIAM Journal on Applied Dynamical Systems 18 (4), 1892-1933, 2019. DOI: https://doi.org/10.1137/18M1225653

and 
Hexagon Invasion Fronts Outside the Homoclinic Snaking Region in the Planar Swift-Hohenberg Equation, DJB Lloyd, 2020, preprint

# Subdirectories:

**IVPs:** Initial value solvers

`swifthohen2DETD_hex_nu_1_6_front_10_defect.m` time evolve a hexagon <10> front that develops a defect
`swifthohen2DETD_hex_invade_patch_fourth_order.m` time evolve a hexagon patch
`swifthohen2DETD_worm_invade_patch_fourth_order.m` time evolve a stripe patch

**BVP_1D_quad_cubic:** `continuation_front_invasion.m`
Continuation of invasion fronts in the quadratic-cubic 1D Swift-Hohenberg equation

**BVP_1D_cubic_quin:** `continuation_front_invasion.m`
Continuation of invasion fronts in the cubic-quintic 1D Swift-Hohenberg equation

**BVP_2D_Stripes:** `continuation_front_invasion.m`
Continuation of almost planar invasion fronts in the cubic-quintic 2D Swift-Hohenberg equation

**BVP_2D_Hex:** `continuation_front_invasion_10_fixed_nonvar.m` and `continuation_front_invasion_11_fixed_nonvar.m` Continuation of hexagon <10> and <11> invasion fronts in the quadratic-cubic 2D Swift-Hohenberg equation

**Hexagon_fold:** `continuation_hexagons_fold_ky_kx_mu_0_32_nu_1_6.m` Continuation of a cellular hexagon fold in k_x and k_y with mu=0.32 and nu=0.9.

**Hexagon_Eckhaus:** `continuation_run_hex_eckhaus.m` Continuation of cellular hexagon Eckhaus boundary in k_x and k_y with mu=0.15 and nu=1.6.




