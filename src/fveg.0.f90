100 continue ! anchor for a GO TO statement that needs to be removed
  !
  NUMIT = 0
  do I = IL1,IL2 ! loop 125
    if (FI(I) > 0. .and. ITER(I) == 1) then
      WZERO(I) = 0.622 * calcEsat(TZERO(I)) / PADRY(I)
      Q0SAT(I) = WZERO(I) / (1.0 + WZERO(I))
      if (IWATER(I) > 0) then
        EVBETA(I) = 1.0
        QZERO(I) = Q0SAT(I)
      else
        !    evaporation coefficient is moss-controlled for peatland--
        if (ipeatland(i) == 0) then
          EVBETA(I) = CEVAP(I)
        else
          ievap(i) = ievapmoss(i)
          evbeta(i) = cevapmoss(i)
        end if
        QZERO(I) = EVBETA(I) * Q0SAT(I) + (1.0 - EVBETA(I)) * QAC(I)
        if (QZERO(I) > QAC(I) .and. IEVAP(I) == 0) then
          EVBETA(I) = 0.0
          QZERO(I) = QAC(I)
        end if
      end if
      !
      TPOTG(I) = TZERO(I) - 8.0 * ZOM(I) * GRAV / SPHAIR
      TVIRTG(I) = TPOTG(I) * (1.0 + 0.61 * QZERO(I))
      if (TVIRTG(I) > TVRTAC(I) + 1.) then
        RAGINV(I) = RAGCO * (TVIRTG(I) - TVRTAC(I)) ** 0.333333
        DRAGIN(I) = 0.333 * RAGCO * (TVIRTG(I) - TVRTAC(I)) ** ( - .667)
      else if (TVIRTG(I) > (TVRTAC(I) + 0.001)) then
        RAGINV(I) = RAGCO * (TVIRTG(I) - TVRTAC(I))
        DRAGIN(I) = RAGCO
      else
        RAGINV(I) = 0.0
        DRAGIN(I) = 0.0
      end if
      !
      QLWOG(I) = SBC * TZERO(I) * TZERO(I) * TZERO(I) * TZERO(I)
      QSENSG(I) = RHOAIR(I) * SPHAIR * RAGINV(I) * &
                  (TPOTG(I) - TAC(I))
      EVAPG (I) = RHOAIR(I) * (QZERO(I) - QAC(I)) * RAGINV(I)
      if (EVAPG(I) > EVPMAX(I)) EVAPG(I) = EVPMAX(I)
      QEVAPG(I) = CPHCHG(I) * EVAPG(I)
      GZERO(I) = GCOEFF(I) * TZERO(I) + GCONST(I)
      RESID(I) = QSWNG(I) + FSVF(I) * QLWIN(I) + (1.0 - FSVF(I)) * &
                 QLWOC(I) - QLWOG(I) - QSENSG(I) - QEVAPG(I) - GZERO(I)
      if (ABS(RESID(I)) < 5.0)                     ITER(I) = 0
      if (ABS(TSTEP(I)) < 1.0E-2)                  ITER(I) = 0
      if (NITER(I) == ITERMX .and. ITER(I) == 1)    ITER(I) = - 1
    end if
  end do ! loop 125
  !
  if (ITCG < 2) then
    !
    !     * OPTION #1: BISECTION ITERATION METHOD.
    !
    do I = IL1,IL2 ! loop 150
      if (FI(I) > 0. .and. ITER(I) == 1) then
        if (NITER(I) == 1) then
          if (RESID(I) > 0.0) then
            TZERO(I) = TZERO(I) + TSTEP(I)
          else
            TZERO(I) = TZERO(I) - TSTEP(I)
          end if
        else
          if ((RESID(I) > 0. .and. TSTEP(I) < 0.) .or. &
              (RESID(I) < 0. .and. TSTEP(I) > 0.)) then
            TSTEP(I) = - TSTEP(I) / 2.0
          end if
          TZERO(I) = TZERO(I) + TSTEP(I)
        end if
        NITER(I) = NITER(I) + 1
        NUMIT = NUMIT + 1
      end if
    end do ! loop 150
    !
  else
    !
    !     * OPTION #2: NEWTON-RAPHSON ITERATION METHOD.
    !
    do I = IL1,IL2 ! loop 175
      if (FI(I) > 0. .and. ITER(I) == 1) then
        DQ0DT = - WZERO(I) * A(I) * (B(I) - TFREZ) / ((TZERO(I) - B(I)) * &
                (1.0 + WZERO(I))) ** 2 * EVBETA(I)
        DRDT0 = - 4.0 * SBC * TZERO(I) ** 3 &
                - GCOEFF(I) - RHOAIR(I) * SPHAIR * &
                (RAGINV(I) + (TPOTG(I) - TAC(I)) * DRAGIN(I)) - &
                CPHCHG(I) * RHOAIR(I) * (DQ0DT * RAGINV(I) &
                + (QZERO(I) - QAC(I)) * DRAGIN(I))
        TSTEP(I) = - RESID(I) / DRDT0
        if (ABS(TSTEP(I)) > 20.0) TSTEP(I) = SIGN(10.0,TSTEP(I))
        TZERO(I) = TZERO(I) + TSTEP(I)
        NITER(I) = NITER(I) + 1
        NUMIT = NUMIT + 1
      end if
    end do ! loop 175
    !
  end if
  !
  if (NUMIT > 0) GO TO 100