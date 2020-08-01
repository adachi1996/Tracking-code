program main
  use const_data
  use map_track
  implicit none
!=====================================================主プログラム
  !磁場マップ使用(map_trackに移動)
  if (select_tracking == 4) then
    call input_map

  !理想磁場使用(mainのまま)
  else
    h           = th_h*(pi/180.0)
    particle(6) = (T**2.0 + 2.0*T*m0)**0.5
    pth         = particle(6)
    !幾何学計算の分岐
    select case (select_B)
      case(1) !Sector型の幾何学計算
        call cal_geo_sec
      case(2) !Rectangular型の幾何学計算
        call cal_geo_rec
      case(3) !工事中
        call cal_geo_rec
      case(4) !FDFトリプレット
        call cal_geo_FDF
    end select
    !計算方法の分岐
    select case (select_tracking)
      case(1) !通常トラッキング(内部サブルーチン)
        call mysub1
      case(2) !しらみつぶし計算(内部サブルーチン)
        call mysub2
      case(3) !アクセプタンス計算(内部サブルーチン)
        call mysub3
      case(5) !安定領域計算(内部サブルーチン)
        call mysub4
    end select
  end if
contains


!==============================================================================================================================
!==============================================================================================================================
!=====================================================[プログラム①:幾何学計算]===================================================
!==============================================================================================================================
!==============================================================================================================================
!================================================================================================プログラムの概要================
!このプログラムは各磁極形状に対応する幾何学計算を行うモジュールである
!いくつか種類があり、const_dataのselect_Bの値によってそれぞれ呼び出されるルーチンが変わる。
!以下に各ルーチンの説明を示す。
!コメントはsector型のみに書いてある。
!各ルーチンで使用されているparticle(8)~(10)はそれぞれBr,Bth,Bzに対応している。
!select_B == 1 : cal_geo_sec
!  初めに設定するパラメータはthF(1/2セルの曲げ角),beta(1/2セルの見込み角),beF・beD(FとDの見込み角),r0(初期入射位置r[m]),セル数N
!  VFFAの計算は主にこのルーチンを使用している。
!select_B == 2 : cal_geo_lec
!  初めに設定するパラメータはbeta(1/2セルの見込み角),beF(Fの見込み角),r0(初期入射位置r[m]),セル数N,L1・L2(磁極の寸法)
!  座標系はlectangular磁石中心に中心を持ち、進行方向をy,外側はxとしている(はず)
!select_B == 3 : cal_geo_FDF
!  工事中('20/3/16現在)
!select_B == 4 : cal_geo_semiFDF
!  何となく作ってみたやつ。少し変わった定義をしていた気がする(足立ログノートNo.2参照)
!  なんとかしてストレートセクションを確保しようとしたパターン。
!  全く使われていないし、今後使うことはないだろう。
!=====================================================[サブルーチン①:Sector型]=====================================================
  subroutine cal_geo_sec  ! 外部サブルーチン
    double precision :: beS
    beS = beta - (beF + beD)
    thD = thF - beta
    roF = r0/(1-cos(thF)+(sin(thF)/tan(beF)))
    r1  = roF*sin(thF)/sin(beF)
    r2  = ((r0 - roF)**2.0-r1**2.-roF**2.0)/(2.0*((r0 - roF)*cos(beF+beS)-r1*cos(beS)))
    roD = r2*sin(beD)/sin(thD)
    r3  = r2*cos(beD) - roD*(1.0 - cos(thD))
    BF = pth/(c*roF*1.0d-6)
    BD = pth/(c*roD*1.0d-6)
    !LF = roF*sin(thF)
    !LD = roD*sin(thD)
    RF = roF*(thF - sin(thF)*cos(thF))/(2.0*sin(thF)) + r1*cos(beF)
    RD = r2*cos(bD) - roD*(thD - sin(thD)*cos(thD))/(2.0*sin(thD))
    if (select_tracking /= 5) then
      print *, "Result of Geometric calculation"
      print "(a,e16.5,e16.5)", "   th0, r0 =",0.0    *180.0/pi, r0
      print "(a,e16.5,e16.5)", "   th1, r1 =",beF    *180.0/pi, r1
      print "(a,e16.5,e16.5)", "   th2, r2 =",beF+beS*180.0/pi, r2
      print "(a,e16.5,e16.5)", "   th3, r3 =",beta   *180.0/pi, r3
      print "(a,e16.5,e16.5)", "   BF , BD =",BF, BD
      print "(a,e16.5)"      , "   FD      =",BF/BD
      print *, ""
    end if
  end subroutine
!=====================================================[サブルーチン②:Lectanguler型]=====================================================
  subroutine cal_geo_rec  ! 外部サブルーチン
    double precision :: tmp1,tmp2,tmp3
    beF_rec  = beF
    thF_rec  = acos((1.0-(r0/L2 - 1.0/tan(beF_rec))**2.)/(1.0+(r0/L2 - 1.0/tan(beF_rec))**2.))
    thD_rec  = thF_rec - beta
    roF      = L2/sin(thF_rec)
    roD      = roF*sin(thF_rec)/sin(thD_rec)
    r1       = L2/sin(beF_rec)
    tmp1     = (r0-roF)**2. - r1**2. - roF**2.
    tmp2     = 2.*L2*((r0-roF)*sin(beta) - r1*sin(beta-beF_rec))
    tmp3     = 2.*   ((r0-roF)*cos(beta) - r1*cos(beta-beF_rec))
    beD_rec  = atan(L2*tmp3/(tmp1-tmp2))
    beS_rec  = beta - beF_rec - beD_rec
    r2       = L2/sin(beD_rec)
    r3       = r2*cos(beD_rec)-roD*(1.0 - cos(thD_rec))
    BF       = pth/(1.0d-6*c*roF)
    BD       = pth/(1.0d-6*c*roD)
    RF       = r1*cos(beF_rec) + 0.5*(L2/sin(thF_rec))*(thF_rec/sin(thF_rec) - cos(thF_rec))
    RD       = r2*cos(beD_rec) - 0.5*(L2/sin(thD_rec))*(thD_rec/sin(thD_rec) - cos(thD_rec))
    if (select_tracking /= 5) then
      write(*,*) beF_rec*180.0/pi,thF_rec*180.0/pi,BF,BD,BF/BD,RF,RD
      write(*,*) r0,r1,r2,r3
    end if
  end subroutine
!=====================================================[サブルーチン③:.FDFトリプレット型]=====================================================
  subroutine cal_geo_FDF  ! 外部サブルーチン
    double precision :: beS
    beta = pi/N
    beS  = beta - (beF + beD)
    thD = thF - beta
    roD = r0/(cos(thD)-1.0+sin(thD)/tan(beD))
    r1  = roD*sin(thD)/sin(beD)
    roF = r1*(sin(beF+beS)-tan(beS)*cos(beF+beS))/(tan(beS)*(1.0-cos(thF))+sin(thF))
    r2  = (r1*sin(beF+beS)-roF*sin(thF))/sin(beS)
    r3  = r2*cos(beS)
    BF  = pth/(1.0d-6*c*roF)
    BD  = pth/(1.0d-6*c*roD)
    Ldr = r2*sin(beS)
    if (select_tracking /= 5) then
      write(*,*) beF*180.0/pi,thD*180.0/pi,BF,BD,BF/BD
      write(*,*) r0,r1,r2,r3
    end if
  end subroutine
!=====================================================[サブルーチン④:エセFDFトリプレット型]=====================================================
  subroutine cal_geo_semiFDF  ! 外部サブルーチン
    double precision :: A,B
    double precision :: L_1,L_2,L_3,L_4
    thD = 2.0*thF - beta
    roF = L2/tan(thF) - L1
    roD = 2.0*L2/tan(thD) - L1
    r1  = ((roD*sin(thD))**2.0 + (roD*(1.0-cos(thD))+r0)**2.0)**0.5
    beD_FDF = asin((roD*sin(thD))/r1)
    roDp= roD + Ls/tan(thD)
    r0p = r0+roD*(1.0-1.0/cos(thD))
    L_1 = (roD*tan(thD))**2.0 + 2.0*Ls*roD*tan(thD) + r1**2.0 - r0p**2.0
    A   = L_1*sin(beD_FDF) - 2.0*roDp*sin(thD)*(r1-r0p*cos(beD_FDF))
    B   = 2.0*roD*r0p*sin(thD)*sin(beD_FDF) - L_1*cos(beD_FDF)
    beS_FDF = atan(A/B)
    r2  = roDp*sin(thD)/sin(beD_FDF+beS_FDF)
    L_2 = roF + (roDp+roF)*sin(thD)/sin(beta)
    L_3 = r0 + roD + Ls/sin(thD) - (roDp + roF)*(cos(thD) + sin(thD)/tan(beta))
    L_4 = L_2**2.0 - L_3**2.0 - 2.0*(roF**2.0)*(1.0-cos(2.0*thF)) + r2**2.0
    A   = 2.0*L_2*sin(beta)*(r2-L_3*cos(beD_FDF+beS_FDF)) - L_4*sin(beD_FDF+beS_FDF)
    B   = 2.0*L_2*L_3*sin(beta)*sin(beD_FDF+beS_FDF) - L_4*cos(beD_FDF+beS_FDF)
    beF_FDF = atan(-A/B)
    r3  = L_2*sin(beta)/sin(beD_FDF+beS_FDF+beF_FDF)
    Ldr = r3*sin(beta-beF_FDF-beD_FDF-beS_FDF)
    BF  = pth/(c*roF*1.0d-6)
    BD  = pth/(c*roD*1.0d-6)
    Lpp = L_2
    Lppp= L_3
    if (select_tracking /= 5) then
      write(*,*) beF_FDF*180.0/pi,thD*180.0/pi,BF,BD,BF/BD
      write(*,*) r0,r1,r2,r3
    end if
  end subroutine

!==============================================================================================================================
!==============================================================================================================================
!=====================================================[プログラム②:軌道計算]=====================================================
!==============================================================================================================================
!==============================================================================================================================
!=====================================================[サブルーチン①:データ保存と軌道計算コード呼び出し]=====================================================通常トラッキング
  subroutine mysub1
    integer :: count_th,count_dth,count_one_cell                              !角度(float)に変わるint型の変数
    double precision :: deruta_P , deruta_Pr, deruta_Pth , RF_T, RF_P, RF_th  !加速関係
    double precision :: temp_deg                                              !角度用の入れ物

    !ファイル作成
    open(17,file=save_name, status='replace')
    write (17,*) 't[s]',',','th[deg]',',','r[m]',',','z[m]',',',' &
                 Pr[MeV/c]',',','Pth[MeV/c]',',','Pz[MeV/c]',',','Br[T]',',','Bth[T]',',','Bz[T]'   !保存データの名前書き込み

    !各変数の初期化
    RF_T = T                                           !加速を考慮したときに使用する運動エネルギーの仮入れ物みたいなもの
    count_th  = 0                                      !角度の代わりに使用するカウント数の初期化
    count_dth = nint(dth/th_h)                         !保存する角度の整数化
    count_one_cell = nint((beta*2.0*180.0/pi)/th_h)    !1cellの角度の整数化

    !入射条件の表示
    print *, "Initial value"
    print *, "  th0 =",particle(2),"[deg.]"
    print *, "  r0  =",particle(3),"[m]"
    print *, "  z0  =",particle(4),"[m]"

    !任意の周回数回るまで計算するループ
    do while (particle(2) <= max_deg)
      if (mod(count_th , count_dth) == 0) then                                                 !任意の角度の整数倍の時にtrue
        write (17,'(e16.8 , a , e16.8 , a , e16.8 , a , e16.8 , a , e16.8 , a , &
                     e16.8 , a , e16.8 , a , e16.8 , a , e16.8 , a , e16.8)') &
                     particle(1),',',particle(2),',',particle(3),',', &
                     particle(4),',',particle(5),',',particle(6),',', &
                     particle(7),',',particle(8),',',particle(9),',',particle(10)                   !データを逐次保存

        !1周したことを知らせるためのif文
        if (mod(count_th,nint(360/th_h)) == 0) then
          print *, nint(particle(2)/360),'Turns'
        end if
      end if

      !理想磁場計算用の角度(0~1cellの角度を繰り返す角度)の計算
      theta = mod(count_th,count_one_cell)*th_h*pi/180.0d0

      !加速関連
      !加速なしの場合でも一応if文の中には入る(なんとなくそうしている)
      if (mod(count_th,nint(360/m_th_h))*m_th_h == th_RF) then
        !加速前の運動量の表示
        print '(a,e16.8,e16.8,e16.8,e16.8)','P(n)   -> ',particle(5),particle(6),particle(7),&
                                                        (particle(5)**2.0 + particle(6)**2.0 + particle(7)**2.0)**0.5
        print '(a,e16.4)','T(n)   -> ',RF_T

        RF_P = (RF_T**2.0 + 2.0*RF_T*m0  -  particle(7)**2.0)**0.5            !加速前の平面上の運動量(スカラ)
        RF_th = ATAN(particle(6)/(-particle(5)))                              !RF_Pの角度
        print '(a,e16.4)','th_1   ->',ATAN(particle(6)/(-particle(5)))*180/pi !加速前の傾きの表示
        RF_T =  RF_T + RF_kV*1d-3
        pth  = (RF_T**2.0 + 2.0*RF_T*m0  -  particle(7)**2.0)**0.5            !加速後の平面上の運動量(スカラ)
        deruta_Pr   = -(pth - RF_P)*cos(RF_th)
        deruta_Pth  =  (pth - RF_P)*sin(RF_th)
        particle(5) = -RF_P*cos(RF_th) + deruta_Pr
        particle(6) =  RF_P*sin(RF_th) + deruta_Pth
        print '(a,e16.4)','th_2   ->',ATAN(particle(6)/(-particle(5)))*180/pi !加速後の傾きの表示
        !加速後の運動量の表示
        print '(a,e16.8,e16.8,e16.8,e16.8)','P(n+1) -> ',particle(5),particle(6),particle(7),&
                                                        (particle(5)**2.0 + particle(6)**2.0 + particle(7)**2.0)**0.5
        print '(a,e16.4)','T(n+1) -> ', -m0+(m0**2.0+particle(5)**2.0 + particle(6)**2.0 + particle(7)**2.0)**0.5
        print *,'------------------------------------------------------'
      end if

      !RKを呼ぶための外部サブルーチンを呼ぶ
      call track
      count_th = count_th + 1      !カウント数の更新
      particle(2) = th_h*count_th  !角度の更新
    end do
    close(17)                          !ファイルを閉じる

    !終了時の表示
    print *, "Last value"
    print *, "  th1 =",th_h*(count_th-1),"[deg.]" !「-1」しているのはDo whileの最後で角度が更新されているため
    print *, "  r1  =",particle(3),"[m]"
    print *, "  z1  =",particle(4),"[m]"
  end subroutine
!=====================================================[サブルーチン②:しらみつぶしプログラム]=====================================================
!角度を考慮できていない初期に作成した閉軌道導出ルーチン
!現在はほとんど使用していない
  subroutine mysub2
    double precision :: temp_r1 = r0 , temp_r2 , ini_r                                              !ｒの情報の入れ物
    double precision :: temp_z1 = z0 , temp_z2 , ini_z                                              !ｚの情報の入れ物
    double precision :: temp_min = 10000.0                                                          !入り口と出口の差の入れ物
    double precision :: temp_deg                                                                    !角度用の入れ物
    integer :: i,j,k                                                                                !整数変数

    do i = 1 , 10                                                                                   !小数点以下iまで計算するループ
      do j = -10 , 10                                                                               !ｒを変化させるためのループ
        ini_r = temp_r1 + dble(j)*dk                                                                !ｒの初期値を更新
        do  k = -10 , 10                                                                            !ｚを変化させるためのループ
          ini_z = temp_z2 + dble(k)*dk                                                              !ｚの初期値を更新
          particle = (/0.0d0, th0, ini_r, ini_z, 0.0d0, pth, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)           !粒子の初期情報を更新
          do while (particle(2) <= 720.0d0)                                                         !1周回るまで計算するループ
            temp_deg = nint(particle(2)/th_h)                                                       !角度の整数化
            theta = mod(temp_deg,(beta*2.0*180.0/pi)/th_h)*th_h*pi/180.0                            !1Cell内での角度の更新
            call track                                                                              !RKを呼ぶための外部サブルーチンを呼ぶ
            particle(2) = particle(2) + th_h                                                        !角度の更新
          end do
          if (temp_min > ((ini_r-particle(3))**2. + (ini_z-particle(4))**2.)**0.5 &                 !より差が小さくなった時true
                                                .and. nint(particle(2)) >= 360) then
            temp_min = ((ini_r-particle(3))**2. + (ini_z-particle(4))**2.)**0.5                     !差の更新
            temp_r2   = ini_r                                                                       !rの更新(仮)
            temp_z2   = ini_z                                                                       !zの更新(仮)
            write(*,*) th0,ini_r,ini_z
            write(*,*) particle(2),particle(3),particle(4)
            write(*,*) '==============================================',temp_min
          end if
        end do
      end do
      dk      = dk*0.1                                                                              !刻み幅の更新
      temp_r1 = temp_r2                                                                             !rの更新(確定)
      temp_z1 = temp_z2                                                                             !zの更新(確定)
    write(*,*) i,temp_min,'dr,dz=',temp_r1-r0,temp_z1
    end do
  end subroutine
!=====================================================[サブルーチン③:安定領域(100turn)の計算用]=====================================================
!任意のr-yの範囲を入射条件として計算するルーチン
!現在ほとんど使用していない
  subroutine mysub3
    integer :: i,j,count                                                                            !整数変数
    double precision :: temp_deg                                                                    !角度用の入れ物
    double precision :: temp_dr,temp_dz                                                             !角度用の入れ物
    count = 0
    do i = 0 , 30
      temp_dr = dble(i-15)*0.001
      do j = 0, 100
        temp_dz = dble(j-100)*0.001
        !write(*,*) temp_dr,temp_dz
        particle = (/0.0d0, th0, r0+temp_dr, z0+temp_dz, 0.0d0, pth, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)   !粒子の初期情報を更新
        if (count == 0) then
          do while (particle(2) <= 36000.0d0)                                                       !100周回るまで計算するループ
            temp_deg = nint(particle(2)/th_h)                                                       !角度の整数化
            theta = mod(temp_deg,(beta*2.0*180.0/pi)/th_h)*th_h*pi/180.0                            !1Cell内での角度の更新
            call track                                                                              !RKを呼ぶための外部サブルーチンを呼ぶ
            particle(2) = particle(2) + th_h                                                        !角度の更新
          end do
        end if
        if (particle(2) >= 35999.0d0 .and. particle(2) <= 36001.0d0) then
          count = 1
          write(*,*) temp_dr,temp_dz
        end if
      end do
      count = 0
    end do
  end subroutine
!=====================================================[サブルーチン④:安定領域(100turn)の計算用]=====================================================
!m値-FD比の安定領域を計算するルーチン
  subroutine mysub4
    integer :: i,j                       !整数変数
    double precision :: temp_deg         !角度用の入れ物
    double precision :: temp_m,temp_thF  !角度用の入れ物
    double precision :: temp_beS,pf      !角度用の入れ物
    write(*,*) 'dr_dz =',dr,dz
    write(*,*) 'm',',','thF',',','FD_ratio'
    do i = 1 , 20
      M = dble(i)*0.5
      do j = 0, 20
        thF = dble(j+20)*pi/180.0
        !temp_beS = dble(j)*0.225
        !beF = (pi/180.0)*(11.25 - temp_beS)/2.0
        !beD = (pi/180.0)*(11.25 - temp_beS)/2.0
        !pf  = (11.25 - temp_beS)/11.25
        !幾何学計算の分岐
        select case (select_B)
          case(1) !Sector型の幾何学計算
            call cal_geo_sec
          case(2) !Rectangular型の幾何学計算
            call cal_geo_rec
          case(3) !工事中
            call cal_geo_rec
          case(4) !FDFトリプレット
            call cal_geo_FDF
        end select
        particle = (/0.0d0, th0, r0+dr, z0+dz, 0.0d0, pth, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)   !粒子の初期情報を更新
        do while (particle(2) <= (rev*360.0) .and. check /= 1)                                                       !100周回るまで計算するループ
          temp_deg = nint(particle(2)/th_h)                                                       !角度の整数化
          theta = mod(temp_deg,(beta*2.0*180.0/pi)/th_h)*th_h*pi/180.0                            !1Cell内での角度の更新
          call track                                                                              !RKを呼ぶための外部サブルーチンを呼ぶ
          particle(2) = particle(2) + th_h                                                        !角度の更新
        end do
        check = 0
        if (particle(2) >= (rev*360.0)) then  !100周した時に入るif文、安定解という判断
          write(*,'(e15.8 , a , e15.8 , a , e15.8 , a , e15.8)') M,',',thF*180.0/pi,',',BF/BD,',',particle(2)
        !if (particle(2) >= (rev*360.0 - 1.0) .and. particle(2) <= (rev*360.0 + 1.0)) then  !100周した時に入るif文、安定解という判断
        !  write(*,'(e15.8 , a , e15.8 , a , e15.8)') M,',',thF*180.0/pi,',',BF/BD
        end if
      end do
    end do
  end subroutine
!=====================================================[サブルーチン⑤:4次のルンゲクッタ]=====================================================
  subroutine track  ! 外部サブルーチン
      double precision :: kl(6),ks(6),k1(6),k2(6),k3(6),k4(6),k(6)                        !ルンゲクッタ用変数入れ
      double precision :: P

      kl = (/particle(1),particle(3),particle(4),particle(5),particle(6),particle(7)/)
      ks = kl

      P      = (particle(5)**2.0+particle(6)**2.0+particle(7)**2.0)**0.5                  !運動量(Ps的な感じ)
      r_beta =    P/(m0**2.0+P**2.0)**0.5

      !βが1以上になったら計算を止めるためのif文
      if (r_beta >= 1.0) then
         !print *, r_beta
         check = 1
      end if
      gamma  = 1.0/(1.0-r_beta**2.0)**0.5

      !計算する理想磁場のタイプの選択
      select case (select_B)
        case(1) !Sector型の幾何学計算
          call sub_B_sec
        case(2) !Rectangular型の幾何学計算
          call sub_B_rec
        case(3) !工事中
          call cal_geo_rec
        case(4) !FDFトリプレット
          call sub_B_FDF
        case(5) !双極磁場
          call sub_dipole
      end select

      !ルンゲクッタの計算
      k1 = RK(ks) !一回目
      ks = kl + 0.5*k1 !更新
      k2 = RK(ks) !二回目
      ks = kl + 0.5*k2 !更新
      k3 = RK(ks) !三回目
      ks = kl + 1.0*k3 !更新
      k4 = RK(ks) !四回目

      k = (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
      kl = kl + k

      particle = (/kl(1),particle(2),kl(2),kl(3),kl(4),kl(5),kl(6),particle(8),particle(9),particle(10)/)

      !計算を終了させるために角度を計算しない程度の大きな値にするためのもの
      if (check == 1) then
        particle(2) = 3600000d0
        check = 0
      end if
  end subroutine
!=====================================================[サブルーチン⑥:円筒座標の運動方程式]=====================================================
  function RK(ks) result(klf)                                                                       !外部関数
    double precision :: ks(6),klf(6)
    klf(1) = h*(ks(2)*gamma*m0/ks(5))/c                                                             !t
    klf(2) = h*(ks(2)*ks(4)/ks(5))                                                                  !r
    klf(3) = h*(ks(2)*ks(6)/ks(5))                                                                  !z
    klf(4) = h*((1.0d-6*c*ks(2)/ks(5))*(ks(5)*particle(10) - ks(6)*particle( 9)) + ks(5))           !pr
    klf(5) = h*((1.0d-6*c*ks(2)/ks(5))*(ks(6)*particle( 8) - ks(4)*particle(10)) - ks(4))           !pth
    klf(6) = h*((1.0d-6*c*ks(2)/ks(5))*(ks(4)*particle( 9) - ks(5)*particle( 8)))                   !pz
  end function


!==============================================================================================================================
!==============================================================================================================================
!=====================================================[プログラム③:理想磁場計算部]===============================================
!==============================================================================================================================
!==============================================================================================================================
  !================================================================================================プログラムの概要==============
  !このプログラムは理想磁場を計算する際に使用する。このプログラムはtrack_dataから呼び出される。
  !いくつか種類があり、const_dataのselect_Bの値によってそれぞれ呼び出されるルーチンが変わる。
  !以下に各ルーチンの説明を示す。基本的にVFFA用に磁場の式は設定されているため、HFFAで使用したい場合は磁場の式を変える必要がある。
  !コメントはsector型のみに書いてある。
  !各ルーチンで使用されているparticle(8)~(10)はそれぞれBr,Bth,Bzに対応している。
  !select_B == 1 : sector
  !  設計軌道上に存在する座標からみたx,y,z,Bx,By,Bzを幾何学計算から求めている。
  !  VFFAの計算は主にこのルーチンを使用している。
  !select_B == 2 : Lectangular
  !  sectorと同様に幾何学計算から磁極のある範囲を計算し、r,thの値から範囲内外の判定をしている。
  !  座標系はlectangular磁石中心に中心を持ち、進行方向をy,外側はxとしている(はず)
  !  漏れ磁場を考慮できていない('20/3/16現在)ため、垂直方向のベータトロン振動の評価には不向き。
  !  最近はほとんど利用していない。2018年くらいまで使用していた(足立ログノートNo.1のどっかにある)
  !select_B == 3 : radial sector
  !  最もシンプルなパターン
  !  単純にr,thで磁石の範囲内外の判定をしている。
  !  漏れ磁場を考慮できていない('20/3/16現在)ため、垂直方向のベータトロン振動の評価には不向き。
  !  最近はほとんど利用していない。2018年くらいまで使用していた(足立ログノートNo.1のどっかにある)
  !select_B == 4 : FDF triplet
  !  何となく作ってみたやつ。少し変わった定義をしていた気がする(足立ログノートNo.2参照)
  !  なんとかしてストレートセクションを確保しようとしたパターン。
  !  全く使われていないし、今後使うことはないだろう。
!=====================================================[サブルーチン①:Sector型]=====================================================
  subroutine sub_B_sec  ! 外部サブルーチン
    double precision :: LF,LD       !各セクター中心からの距離を入れる変数
    double precision :: temp1,temp2 !幾何学計算用仮入れ物
    double precision :: faiF,faiD   !セクター中心からの角度を入れる
    double precision :: r,z,x       !トラッキングパラメータ

    r = particle(3)
    z = particle(4)

    !rやzがあまりにも設計軌道から外れていた場合は計算をストップしたかったif文
    !値に科学的な根拠はなく、経験？というか何となくで決めた
    if (abs(r) > (r0 + 0.1) .or. abs(r) < (r0 - 0.1) .or. abs(z) > (z0 + 0.01 + dz*2.0)) then
      check = 1
    end if

    !1セル前半(1/2)の計算
    if (theta <= beta) then
      temp1 = r**2.0 - 2.0*r*(r0-roF)*cos(theta)+(r0-roF)**2.0        !F側のセクター中心からの距離を計算
      LF   = temp1**0.5                                               !temp1の値が知りたかったためあえてこのようにしている
      faiF = asin((r*sin(theta))/LF)                                  !セクター中心からの角度を計算

      temp2 = r**2.0 - 2.0*r*(r3+roD)*cos(beta - theta)+(r3+roD)**2.0 !F側のセクター中心からの距離を計算(以下同文)
      LD   = temp2**0.5
      faiD = asin((r*sin(beta - theta))/LD)

      !数学的におかしな数値になった場合は計算をストップする(temp1やtemp2はここで使用)
      if (temp1 <= 0.0 .or. (r*sin(theta))/LF        > 1.0 .or. &
          temp2 <= 0.0 .or. (r*sin(beta - theta))/LD > 1.0) then
        check = 1
      end if

      if (faiF <= thF .and. faiD > thD) then !F側の計算
        x   = LF - roF !設計軌道からのずれ
        if (abs(M*x) > pi/2.0) then !磁場の定義式の範囲内がどうか判定
          check = 1
        end if
        particle( 8) =  BF*exp(M*z)*sin(M*x)
        particle( 9) =  0.
        particle(10) = -BF*exp(M*z)*cos(M*x)
        !particle( 8) =  BF*M*x
        !particle( 9) =  0.
        !particle(10) = -BF*(1.0+M*z)
      else if (faiD <= thD .and. faiF > thF) then !D側の計算
        x   = roD - LD !設計軌道からのずれ
        if (abs(M*x) > pi/2.0) then !磁場の定義式の範囲内がどうか判定
          check = 1
        end if
        particle( 8) = -BD*exp(M*z)*sin(M*x)
        particle( 9) =  0.
        particle(10) =  BD*exp(M*z)*cos(M*x)
        !particle( 8) = -BD*M*x
        !particle( 9) =  0.
        !particle(10) =  BD*(1.0+M*z)
      else
        particle( 8) =  0.0
        particle( 9) =  0.0
        particle(10) =  0.0
      end if

    else !残り(1/2セル)を計算
      temp1 = r**2.0 - 2.0*r*(r3+roD)*cos(theta - beta)+(r3+roD)**2.0
      LD   = temp1**0.5
      faiD = asin((r*sin(theta - beta))/LD)

      temp2 = r**2.0 - 2.0*r*(r0-roF)*cos(2.0*beta - theta)+(r0-roF)**2.0
      LF   = temp2**0.5
      faiF = asin((r*sin(2.0*beta - theta))/LF)

      if (temp1 <= 0.0 .or. (r*sin(theta - beta))/LD     > 1.0 .or. &
          temp2 <= 0.0 .or. (r*sin(2.0*beta - theta))/LF > 1.0) then
        check = 1
      end if

      if (faiD <= thD .and. faiF > thF) then
        x   =  roD - LD
        if (abs(M*x) > pi/2.0) then
          check = 1
        end if
        particle( 8) = -BD*exp(M*z)*sin(M*x)
        particle( 9) =  0.
        particle(10) =  BD*exp(M*z)*cos(M*x)
        !particle( 8) = -BD*M*x
        !particle( 9) =  0.
        !particle(10) =  BD*(1.0+M*z)

      else if (faiF <= thF .and. faiD > thD) then
        x   =  LF - roF
        if (abs(M*x) > pi/2.0) then
          check = 1
        end if
        particle( 8) =  BF*exp(M*z)*sin(M*x)
        particle( 9) =  0.
        particle(10) = -BF*exp(M*z)*cos(M*x)
        !particle( 8) =  BF*M*x
        !particle( 9) =  0.
        !particle(10) = -BF*(1.0+M*z)
      else
        particle( 8) =  0.0
        particle( 9) =  0.0
        particle(10) =  0.0
      end if
    end if
  end subroutine
!=====================================================[サブルーチン②:Rectanguler型]=====================================================
  subroutine sub_B_rec  ! 外部サブルーチン
    double precision :: LF,LD,x,B
    !print *, theta*180.0/pi,beta*180.0/pi,RF,RD
    if (theta <= beta) then
      LF = particle(3)*sin(theta)
      LD = particle(3)*sin(beta - theta)
      if (LF <= L2 .and. LD <=L2) then
        particle(10) = 1000.0d0
      else if (LF <= L2) then
        x = particle(3)*cos(theta) - RF
        B = BF*exp(M*particle(4))*sin(M*x)
        particle( 8) =  B*cos(theta)
        particle( 9) = -B*sin(theta)
        particle(10) = -BF*exp(M*particle(4))*cos(M*x)
      else if (LD <= L2) then
        x = particle(3)*cos(beta-theta) - RD
        B = -BD*exp(M*particle(4))*sin(M*x)
        particle( 8) =  B*cos(beta-theta)
        particle( 9) =  B*sin(beta-theta)
        particle(10) =  BD*exp(M*particle(4))*cos(M*x)
      else
        particle( 8) = 0.0d0
        particle( 9) = 0.0d0
        particle(10) = 0.0d0
      end if
    else
      LD = particle(3)*sin(theta-beta)
      LF = particle(3)*sin(2.*beta-theta)
      if (LF <= L2 .and. LD <=L2) then
        particle(10) = 1000.0d0
        return
      else if (LD <= L2) then
        x = particle(3)*cos(theta-beta) - RD
        B = -BD*exp(M*particle(4))*sin(M*x)
        particle( 8) =  B*cos(theta-beta)
        particle( 9) = -B*sin(theta-beta)
        particle(10) =  BD*exp(M*particle(4))*cos(M*x)
      else if (LF <= L2) then
        x = particle(3)*cos(2.*beta-theta) - RF
        B = BF*exp(M*particle(4))*sin(M*x)
        particle( 8) =  B*cos(2.*beta-theta)
        particle( 9) = -B*sin(2.*beta-theta)
        particle(10) = -BF*exp(M*particle(4))*cos(M*x)
      else
        particle( 8) = 0.0d0
        particle( 9) = 0.0d0
        particle(10) = 0.0d0
      end if
    end if
  end subroutine
!=====================================================[サブルーチン③:FDFトリプレット型]=====================================================
  subroutine sub_B_FDF  ! 外部サブルーチン
    double precision :: LF,LD,x,faiF,faiD,fai
    double precision :: r,z

    r = particle(3)
    z = particle(4)

    if (theta <= beta) then!1/2セルの前半
      LD   = (r**2.0 - 2.0*r*(r0+roD)*cos(theta)+(r0+roD)**2.0)**0.5
      faiD = asin(r*sin(theta)/LD)
      LF   = (r**2.0 - 2.0*r*(Ldr*sin(beta-theta)-(roF-r3)*cos(beta-theta))+Ldr**2.0+(roF-r3)**2.0)**0.5
      faiF = asin((r*sin(beta-theta)-Ldr)/LF)
      if (faiD <= thD .and. faiD >= 0.0) then!D磁石
        x   = roD - LD
        fai = M*x
        if (abs(fai) > 90.0*pi/180.0) then
          check = 1
        end if
        particle( 8) = -BD*exp(M*particle(4))*sin(fai)
        particle( 9) =  0.0
        particle(10) =  BD*exp(M*particle(4))*cos(fai)

      else if (faiF <= thF .and. r*sin(beta-theta) >= Ldr) then!F磁石
        x   = LF - roF
        fai = M*x
        if (abs(fai) > 90.0*pi/180.0) then
          check = 1
        end if
        particle( 8) =  BF*exp(M*particle(4))*sin(fai)
        particle( 9) =  0.0
        particle(10) = -BF*exp(M*particle(4))*cos(fai)

      else !ストレートセクション
        particle( 8) =  0.0
        particle( 9) =  0.0
        particle(10) =  0.0
      end if
    end if

    if (theta > beta) then!1/2セルの後半
      LF   = (r**2.0 - 2.0*r*(Ldr*sin(theta-beta)-(roF-r3)*cos(theta-beta))+Ldr**2.0+(roF-r3)**2.0)**0.5
      faiF = asin((r*sin(theta-beta)-Ldr)/LF)
      LD   = (r**2.0 - 2.0*r*(r0+roD)*cos(2.0*beta-theta)+(r0+roD)**2.0)**0.5
      faiD = asin(r*sin(2.0*beta-theta)/LD)
      if (faiF <= thF .and. r*sin(theta-beta) >= Ldr) then!D磁石
        x   =  LF - roF
        fai =  M*x
        if (abs(fai) > 90.0*pi/180.0) then
          check = 1
        end if
        particle( 8) =  BF*exp(M*particle(4))*sin(fai)
        particle( 9) =  0.0
        particle(10) = -BF*exp(M*particle(4))*cos(fai)
      else if (faiD <= thD .and. faiD >= 0.0) then !F磁石
        x   =  roD - LD
        fai =  M*x
        if (abs(fai) > 90.0*pi/180.0) then
          check = 1
        end if
        particle( 8) = -BD*exp(M*particle(4))*sin(fai)
        particle( 9) =  0.0
        particle(10) =  BD*exp(M*particle(4))*cos(fai)
      else !ストレートセクション
        particle( 8) =  0.0
        particle( 9) =  0.0
        particle(10) =  0.0
      end if
    end if
  end subroutine
!=====================================================[サブルーチン④:エセFDFトリプレット型]=====================================================
  subroutine sub_B_semiFDF  ! 外部サブルーチン
    double precision :: LD,LF,x,B
    double precision :: faiF,faiD,fai1
    double precision :: r,z

    r = particle(3)
    z = particle(4)

    if (theta <= beta) then !1Cellの半分
      LD   = (r**2.0 - 2.0*(roD+r0)*r*cos(theta) + (roD+r0)**2.0)**0.5
      faiD = asin(r*sin(theta)/LD)

      LF   = ((r*sin(beta - theta)-Ldr)**2.0 + (r*cos(beta - theta)-(Lpp-roF)-Lppp*cos(beta))**2.0)**0.5
      faiF = asin((r*sin(beta - theta)-Ldr)/LF)

      if (faiD >= 0.0 .and. faiD < thD) then
        x = roD - LD
        if (abs(M*x) > pi/2.0) then
          check = 1
        end if
        particle( 8) = -BD*exp(M*z)*sin(M*x)
        particle( 9) =  0.0
        particle(10) =  BD*exp(M*z)*cos(M*x)
      else if (faiF > 0.0 .and. faiF <= 2.0*thF) then
        x = LF - roF
        if (abs(M*x) > pi/2.0) then
          check = 1
        end if
        particle( 8) =  BF*exp(M*z)*sin(M*x)
        particle( 9) =  0.0
        particle(10) = -BF*exp(M*z)*cos(M*x)
      else
        particle( 8) =  0.0
        particle( 9) =  0.0
        particle(10) =  0.0
      end if
    else
      LF   = ((r*sin(theta-beta)-Ldr)**2.0 + (r*cos(theta-beta)-(Lpp-roF)-Lppp*cos(beta))**2.0)**0.5
      faiF = asin((r*sin(theta-beta)-Ldr)/LF)

      LD   = (r**2.0 - 2.0*(roD+r0)*r*cos(2.0*beta - theta) + (roD+r0)**2.0)**0.5
      faiD = asin(r*sin(2.0*beta - theta)/LD)

      if (faiF > 0.0 .and. faiF <= 2.0*thF) then
        x = LF - roF
        if (abs(M*x) > pi/2.0) then
          check = 1
        end if
        particle( 8) =  BF*exp(M*z)*sin(M*x)
        particle( 9) =  0.0
        particle(10) = -BF*exp(M*z)*cos(M*x)
      else if (faiD >= 0.0 .and. faiD < thD) then
        x = roD - LD
        if (abs(M*x) > pi/2.0) then
          check = 1
        end if
        particle( 8) = -BD*exp(M*z)*sin(M*x)
        particle( 9) =  0.0
        particle(10) =  BD*exp(M*z)*cos(M*x)
      else
        particle( 8) =  0.0
        particle( 9) =  0.0
        particle(10) =  0.0
      end if
    end if
    !print *,theta*180.0/pi,faiD*180.0/pi,faiF*180.0/pi,r
  end subroutine
!=====================================================[サブルーチン⑤単純なダイポール:]=====================================================
  subroutine sub_dipole  ! 外部サブルーチン
    particle( 8) =  0.0
    particle( 9) =  0.0
    particle(10) = -6.67128d-05
  end subroutine
end program
