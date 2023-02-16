/*

    Utility: OSCA;
    Module: VQTL;
    Methods: drm, svlm;

    1. drm method:
    C implementation of Deviation Regression Model (DRM).
    The original R code is available at https://github.com/drewmard/DRM.
    See https://doi.org/10.1016/j.ajhg.2020.11.016.

    2. svlm method:
    Implementation of SVLM(Squared residual Value Linear Modeling).
    The R package is available at
    https://cran.r-project.org/src/contrib/Archive/VariABEL/VariABEL_0.9-2.tar.gz.
    See https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-13-4.


    -- Benjamin Fang 20230117
*/

#ifndef OSCA_VQTL_DRM_SVLM_HEAD
#define OSCA_VQTL_DRM_SVLM_HEAD

#if defined OSCA_VQTL_DRM_SVLM_SRC
#define MODULE_VQTL_DRM_EXTERN
#else
#define MODULE_VQTL_DRM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

MODULE_VQTL_DRM_EXTERN int Module_vqtl_drm(int argc, char *argv[]);
MODULE_VQTL_DRM_EXTERN int Module_vqtl_svlm(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif
