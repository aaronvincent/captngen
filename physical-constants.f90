module phys
    !* The physical constants come from the tables found in the
    ! [PDG](https://pdg.lbl.gov/2025/reviews/contents_sports.html#collapseListGroupConstants), where 'PDG 2024' indicates the value
    ! is from the 2024 release of the PDG, and 'exact' inicates that the value is precise to all places. Lacking these, the value is
    ! derived from other constants.
    implicit none

    double precision, parameter :: pi = 3.141592653589793238d0
    !! \( \pi \) (PDG 2024) [\( \text{1} \)]

    double precision, parameter :: GN = 6.67430d-8
    !! Newton's gravitational constant \( G_N \) (PDG 2024) [\( \text{cm}^3 \text{g}^{-1} \text{s}^{-2} \)]

    double precision, parameter :: NAvo = 6.02214076d23
    !! Avogadro's constant \( N_A \) (PDG 2024 exact) [\( \text{mol}^{-1} \)]

    double precision, parameter :: mass_sun = 1.98841d33
    !! Mass of the Sun \( M_\odot \) (PDG 2024) [\( \text{g} \)]

    double precision, parameter :: radius_sun = 6.957d10
    !! Radius of the Sun \( R_\odot \) (PDG 2024 exact) [\( \text{cm} \)]

    double precision, parameter :: GMoverR = GN*mass_sun/radius_sun
    !! \( \frac{ G M_\odot }{ R_\odot } \) [\( \text{cm}^2 \text{s}^{-2} \)]

    double precision, parameter :: c0 = 2.99792458d10
    !! Speed of light \( c \) (PDG 2024 exact) [\( \text{cm} \  \text{s}^{-1} \)]

    double precision, parameter :: mnuc = 0.93827208816d0
    !! Proton mass \( m_P \) (PDG 2024) [\( \text{GeV} \)]

    double precision, parameter :: kB = 1.380649d-16
    !! Boltzmann constant \( k_B \) (PDG 2024 exact) [\( \text{erg} \  \text{K}^{-1} \)]

    double precision, parameter :: electric = 1.602176634d-19
    !! Elementary charge \( e \) (PDG 2024 exact) [\( \text{C} \)]

    double precision, parameter :: gev_erg = 1.d-16/electric
    !! GeV per erg \( \frac{\text{GeV}}{\text{erg}} \) [\( 10^{-9} \text{GeV} = e \times 10^7 \text{erg} \  \text{C}^{-1} \)]

    double precision, parameter :: hbar = 6.62607015d-27/(2.d0*pi) * gev_erg
    !! Reduced Planck's constant \( \hbar \) (PDG 2024) [\( \text{GeV} \  \text{s} \)]

end module phys
