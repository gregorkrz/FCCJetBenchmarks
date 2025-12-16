Main:numberOfEvents = 50000
Main:timesAllowErrors = 5

Init:showChangedSettings = on
Init:showChangedParticleData = off
Next:numberCount = 100
Next:numberShowInfo = 1
Next:numberShowProcess = 1
Next:numberShowEvent = 0

! --- e+e- beams
Beams:idA = 11
Beams:idB = -11
Beams:eCM = 240


! --- Hard process: e+ e- -> Z H (Higgsstrahlung)
HiggsSM:ffbar2HZ = on

! --- Showering
PartonLevel:ISR = on
PartonLevel:FSR = on

! --- Higgs mass
25:m0 = 125.0

! --- Force Z hadronic decays only: Z -> qqbar (q = d,u,s,c,b)
23:onMode = off
23:onIfAny = 1 2 3 4 5

! --- Force Higgs -> W+ W- (off-shell allowed at mH=125)
25:onMode = off
25:onIfMatch = 24 -24

! --- Force both W's to hadronic decays only: W -> q q'  (exclude leptons)
24:onMode = off
24:onIfAny = 1 2 3 4 5

! --- Switch on QED FSR from charged particles (e+/e- and quarks)
