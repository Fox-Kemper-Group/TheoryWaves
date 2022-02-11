!BOP

! !IROUTINE: cvmix_kpp_compute_unresolved_shear
! !INTERFACE:

  function cvmix_kpp_compute_unresolved_shear(zt_cntr, ws_cntr, N_iface,      &
                                            Nsqr_iface, EFactor,              &
                                            LaSL, bfsfc, ustar,               &
                                            CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the square of the unresolved shear ($V_t^2$ in Eq. (23) of LMD94)
!  at cell centers. Note that you must provide either the buoyancy frequency
!  or its square at cell interfaces, this routine by default will use the
!  lower cell interface value as the cell center, but you can instead take
!  an average of the top and bottom interface values by setting
!  lavg\_N\_or\_Nsqr = .true. in cvmix\_kpp\_init(). If you pass in Nsqr then
!  negative values are assumed to be zero (default POP behavior).
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    ! zt_cntr: height at center of cell (units: m)
    ! ws_cntr: w_s (turbulent scale factor) at center of cell (units: m/s)
    real(cvmix_r8), dimension(:), intent(in) :: zt_cntr,  ws_cntr
    ! N_iface: buoyancy frequency at cell interfaces (units: 1/s)
    ! Nsqr_iface: squared buoyancy frequency at cell interfaces (units: 1/s^2)
    ! note that you must provide exactly one of these two inputs!
    real(cvmix_r8), dimension(size(zt_cntr)+1), intent(in), optional ::       &
                                                    N_iface, Nsqr_iface
    ! EFactor: Langmuir enhancement factor (units: none)
    ! LaSL: surface layer averaged Langmuir number (units: none)
    ! bfsfc: surface buoyancy flux (units: m^2/s^3)
    ! ustar: friction velocity (units: m/s)
    real(cvmix_r8), intent(in), optional :: EFactor, LaSL, bfsfc, ustar
    type(cvmix_kpp_params_type),  intent(in), optional, target ::             &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(size(zt_cntr)) ::                               &
                             cvmix_kpp_compute_unresolved_shear

!EOP
!BOC

    ! Local variables
    integer :: kt, nlev
    real(cvmix_r8) :: Cv, Vtc
    ! N_cntr: buoyancy frequency at cell centers, derived from either N_iface
    !        or Nsqr_iface (units: 1/s)
    real(cvmix_r8), dimension(size(zt_cntr)) :: N_cntr
    ! c_CT, c_ST, c_LT, p_LT: parameters of Langmuir-enhanced entrainment
    !                         in Li and Fox-Kemper, 2017, JPO
    real(cvmix_r8) :: c_CT, c_ST, c_LT, p_LT
    ! RWHGK_ENTR_COEF, RWHGK_ENTR_EXP: parameters of Langmuir-enhanced
    !                         entrainment in Reichl et al., 2016, JPO
    real(cvmix_r8) :: RWHGK_ENTR_COEF, RWHGK_ENTR_EXP
    ! Vt2_Enhancement: enhancement factor for unresolved shear
    real(cvmix_r8) :: Vt2_Enhancement
    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_in

    nlev = size(zt_cntr)
    if (size(ws_cntr).ne.nlev) then
      print*, "ERROR: zt_cntr and ws_cntr must be same size"
      stop 1
    end if

    if (present(N_iface).and.present(Nsqr_iface)) then
      print*, "ERROR: you must provide N_iface OR Nsqr_iface, can not send",  &
              "both!"
      stop 1
    end if

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    if (present(N_iface)) then
      if (size(N_iface).ne.(nlev+1)) then
        print*, "ERROR: N_iface must have one more element than zt_cntr"
        stop 1
      end if
      do kt=1,nlev
        N_cntr(kt) = N_iface(kt+1)
      end do
    else
      if (present(Nsqr_iface)) then
        if (size(Nsqr_iface).ne.(nlev+1)) then
          print*, "ERROR: Nsqr_iface must have one more element than zt_cntr"
          stop 1
        end if
        do kt=1,nlev
          N_cntr(kt)=sqrt(max(Nsqr_iface(kt+1),cvmix_zero))
        end do
      else
        print*, "ERROR: you must provide N_iface OR Nsqr_iface"
        stop 1
      end if
    end if

    ! options for Langmuir enhanced entrainment
    select case (CVmix_kpp_params_in%Langmuir_Entrainment_Opt)

      case (LANGMUIR_ENTRAINMENT_LWF16)
        if (.not.(present(EFactor) )) then
           print*, "ERROR: you must pass in EFactor if ",&
                "Langmuir_entrainment_str .eq. 'LWF16'!"
           stop 1
        end if
        Vt2_Enhancement = EFactor

        ! From LMD 94, Vtc = sqrt(-beta_T/(c_s*eps))/kappa^2
        Vtc = sqrt(0.2_cvmix_r8/(cvmix_get_kpp_real('c_s', CVmix_kpp_params_in) * &
              cvmix_get_kpp_real('surf_layer_ext', CVmix_kpp_params_in))) / &
              (cvmix_get_kpp_real('vonkarman', CVmix_kpp_params_in)**2)

        do kt=1,nlev
          if (CVmix_kpp_params_in%lscalar_Cv) then
            Cv = cvmix_get_kpp_real('Cv', CVmix_kpp_params_in)
          else
            ! Cv computation comes from Danabasoglu et al., 2006
            if (N_cntr(kt).lt.0.002_cvmix_r8) then
              Cv = 2.1_cvmix_r8-real(200,cvmix_r8)*N_cntr(kt)
            else
              Cv = 1.7_cvmix_r8
            end if
          end if

          cvmix_kpp_compute_unresolved_shear(kt) = -Cv*Vtc*zt_cntr(kt)*          &
                                N_cntr(kt)*ws_cntr(kt)/                          &
                                CVmix_kpp_params_in%Ri_crit * Vt2_Enhancement
          if (cvmix_kpp_compute_unresolved_shear(kt).lt.                         &
              CVmix_kpp_params_in%minVtsqr) then
            cvmix_kpp_compute_unresolved_shear(kt) = CVmix_kpp_params_in%minVtsqr
          end if
        end do

      case (LANGMUIR_ENTRAINMENT_LF17)

        if (.not.(present(LaSL) .and. present(bfsfc) .and. present(ustar))) then
          print*, "ERROR: you must pass in LaSL, bfsfc and ustar if ",&
                "Langmuir_entrainment_str == 'LF17'!"
          stop 1
        end if
        ! only apply Langmuir enhanced entrainment under unstable condition
        if (bfsfc<cvmix_zero) then
          ! (26) of Li and Fox-Kemper, 2017, JPO
          c_CT =  cvmix_get_kpp_real('c_CT', CVmix_kpp_params_in)
          c_ST =  cvmix_get_kpp_real('c_ST', CVmix_kpp_params_in)
          c_LT =  cvmix_get_kpp_real('c_LT', CVmix_kpp_params_in)
          p_LT =  cvmix_get_kpp_real('p_LT', CVmix_kpp_params_in)
          do kt=1,nlev
            if (CVmix_kpp_params_in%lscalar_Cv) then
              Cv = cvmix_get_kpp_real('Cv', CVmix_kpp_params_in)
            else
            ! Cv computation comes from Danabasoglu et al., 2006
              if (N_cntr(kt).lt.0.002_cvmix_r8) then
                Cv = 2.1_cvmix_r8-real(200,cvmix_r8)*N_cntr(kt)
              else
                Cv = 1.7_cvmix_r8
              end if
            end if
            Vtc = sqrt((c_CT*bfsfc*zt_cntr(kt) + c_ST*ustar**3 +                 &
                     c_LT*ustar**3*LaSL**(-1.*p_LT))/ws_cntr(kt))
            cvmix_kpp_compute_unresolved_shear(kt) = -Cv*Vtc*zt_cntr(kt)*        &
                                   N_cntr(kt)/CVmix_kpp_params_in%Ri_crit
            if (cvmix_kpp_compute_unresolved_shear(kt).lt.                       &
                CVmix_kpp_params_in%minVtsqr) then
              cvmix_kpp_compute_unresolved_shear(kt) = CVmix_kpp_params_in%minVtsqr
            end if
          end do
        else
          ! From LMD 94, Vtc = sqrt(-beta_T/(c_s*eps))/kappa^2
          Vtc = sqrt(0.2_cvmix_r8/(cvmix_get_kpp_real('c_s', CVmix_kpp_params_in) * &
                cvmix_get_kpp_real('surf_layer_ext', CVmix_kpp_params_in))) / &
                (cvmix_get_kpp_real('vonkarman', CVmix_kpp_params_in)**2)

          do kt=1,nlev
            if (CVmix_kpp_params_in%lscalar_Cv) then
              Cv = cvmix_get_kpp_real('Cv', CVmix_kpp_params_in)
            else
              ! Cv computation comes from Danabasoglu et al., 2006
              if (N_cntr(kt).lt.0.002_cvmix_r8) then
                Cv = 2.1_cvmix_r8-real(200,cvmix_r8)*N_cntr(kt)
              else
                Cv = 1.7_cvmix_r8
              end if
            end if

            cvmix_kpp_compute_unresolved_shear(kt) = -Cv*Vtc*zt_cntr(kt) *       &
                              N_cntr(kt)*ws_cntr(kt)/CVmix_kpp_params_in%Ri_crit
            if (cvmix_kpp_compute_unresolved_shear(kt).lt.                       &
                CVmix_kpp_params_in%minVtsqr) then
              cvmix_kpp_compute_unresolved_shear(kt) = CVmix_kpp_params_in%minVtsqr
            end if
          end do
        end if

      case (LANGMUIR_ENTRAINMENT_RWHGK16)

        if (.not.(present(LaSL) )) then
           print*, "ERROR: you must pass in LaSL if ",&
                "Langmuir_entrainment_str == 'RWHGK16'!"
           stop 1
        end if
        RWHGK_ENTR_COEF =  cvmix_get_kpp_real('RWHGK_ENTR_COEF', &
             CVmix_kpp_params_in)
        RWHGK_ENTR_EXP =  cvmix_get_kpp_real('RWHGK_ENTR_EXP', &
             CVmix_kpp_params_in)
        Vt2_Enhancement = cvmix_one + RWHGK_ENTR_COEF * LASL**RWHGK_ENTR_EXP

        ! From LMD 94, Vtc = sqrt(-beta_T/(c_s*eps))/kappa^2
        Vtc = sqrt(0.2_cvmix_r8/(cvmix_get_kpp_real('c_s', CVmix_kpp_params_in) * &
              cvmix_get_kpp_real('surf_layer_ext', CVmix_kpp_params_in))) / &
              (cvmix_get_kpp_real('vonkarman', CVmix_kpp_params_in)**2)

        do kt=1,nlev
          if (CVmix_kpp_params_in%lscalar_Cv) then
            Cv = cvmix_get_kpp_real('Cv', CVmix_kpp_params_in)
          else
            ! Cv computation comes from Danabasoglu et al., 2006
            if (N_cntr(kt).lt.0.002_cvmix_r8) then
              Cv = 2.1_cvmix_r8-real(200,cvmix_r8)*N_cntr(kt)
            else
              Cv = 1.7_cvmix_r8
            end if
          end if

          cvmix_kpp_compute_unresolved_shear(kt) = -Cv*Vtc*zt_cntr(kt)*          &
                                N_cntr(kt)*ws_cntr(kt)/                          &
                                CVmix_kpp_params_in%Ri_crit * Vt2_Enhancement
          if (cvmix_kpp_compute_unresolved_shear(kt).lt.                         &
              CVmix_kpp_params_in%minVtsqr) then
            cvmix_kpp_compute_unresolved_shear(kt) = CVmix_kpp_params_in%minVtsqr
          end if
        end do

      case DEFAULT

        ! From LMD 94, Vtc = sqrt(-beta_T/(c_s*eps))/kappa^2
        Vtc = sqrt(0.2_cvmix_r8/(cvmix_get_kpp_real('c_s', CVmix_kpp_params_in) * &
              cvmix_get_kpp_real('surf_layer_ext', CVmix_kpp_params_in))) / &
              (cvmix_get_kpp_real('vonkarman', CVmix_kpp_params_in)**2)

        do kt=1,nlev
          if (CVmix_kpp_params_in%lscalar_Cv) then
            Cv = cvmix_get_kpp_real('Cv', CVmix_kpp_params_in)
          else
            ! Cv computation comes from Danabasoglu et al., 2006
            if (N_cntr(kt).lt.0.002_cvmix_r8) then
              Cv = 2.1_cvmix_r8-real(200,cvmix_r8)*N_cntr(kt)
            else
              Cv = 1.7_cvmix_r8
            end if
          end if

          cvmix_kpp_compute_unresolved_shear(kt) = -Cv*Vtc*zt_cntr(kt) *       &
                            N_cntr(kt)*ws_cntr(kt)/CVmix_kpp_params_in%Ri_crit
          if (cvmix_kpp_compute_unresolved_shear(kt).lt.                       &
              CVmix_kpp_params_in%minVtsqr) then
            cvmix_kpp_compute_unresolved_shear(kt) = CVmix_kpp_params_in%minVtsqr
          end if
        end do

    end select

!EOC

  end function cvmix_kpp_compute_unresolved_shear
