C> \ingroup nwxc
C> @{
C>
C> \file nwxc_block.F
C> Functional initialization
C>
C> @}
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief One time initialization of the density functional library
C> data structures
C>
C> Initialize the density functional library data structures ensuring
C> that initially there are no functional terms and all weights are
C> zero.
C>
      block data nwxc_block
      implicit none
#include "nwxcP.fh"
#include "nwxc_param.fh"
      data nwxc_depends_on_rho   / .false. /
      data nwxc_depends_on_gamma / .false. /
      data nwxc_depends_on_tau   / .false. /
      data nwxc_cam_srhfex       / .false. /
      data nwxc_on               / .false. /
      data nwxc_num_f      / 0 /
      data nwxc_usr_num    / 0 /
      data nwxc_wght_hf    / 0.0d0 /
      data nwxc_wght_mp2   / 0.0d0 /
      data nwxc_cam_alpha  / 0.0d0 /
      data nwxc_cam_beta   / 0.0d0 /
      data nwxc_cam_gamma  / 0.0d0 /
      data nwxc_rho_tol    / 1.0d-18 /
      data nwxc_disp_id    / NWXC_DISP_NONE /
      data nwxc_disp_s6    / 0.0d0 /
      data nwxc_disp_s8    / 0.0d0 /
      data nwxc_disp_alpha / 0.0d0 /
      end
C> @}
C> \ingroup nwxc_api
C> @{
C>
C> \brief Initialize the density functional library data structures
C>
C> Explicitly initialize the density functional library data structures.
C> This initialization is designed to match that of the `nwxc_block`
C> block data code. The only difference is that this routine can be
C> invoked explicitly at a moment of the application programmers choice.
C>
      subroutine nwxc_init
      implicit none
#include "nwxcP.fh"
#include "nwxc_param.fh"
      nwxc_depends_on_rho   = .false.
      nwxc_depends_on_gamma = .false.
      nwxc_depends_on_tau   = .false.
      nwxc_cam_srhfex       = .false.
      nwxc_on               = .false. 
      nwxc_num_f      = 0 
      nwxc_usr_num    = 0 
      nwxc_wght_hf    = 0.0d0 
      nwxc_wght_mp2   = 0.0d0
      nwxc_cam_alpha  = 0.0d0
      nwxc_cam_beta   = 0.0d0
      nwxc_cam_gamma  = 0.0d0
      nwxc_rho_tol    = 1.0d-18
      nwxc_disp_id    = NWXC_DISP_NONE
      nwxc_disp_s6    = 0.0d0
      nwxc_disp_s8    = 0.0d0
      nwxc_disp_alpha = 0.0d0
      end
C> @}
