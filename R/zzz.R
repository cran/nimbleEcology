# This file registers all distributions when the package is loaded.
.onAttach <- function(libname, pkgname) {

  packageStartupMessage("Loading nimbleEcology. Registering multiple variants of the following distributions:\n ",
                        "dOcc", ", dDynOcc", ", dCJS", ", dHMM", ", dDHMM", ", dNmixture.\n")

# Register the distributions explicitly for two reasons:
# 1. Avoid message to user about automatic registrations upon first use in a nimbleModel
# 2. Establish default len = 0 via reparameterization mechanism.
  suppressMessages({
  registerDistributions(list(
    dCJS_ss = list(
      BUGSdist = "dCJS_ss(probSurvive, probCapture, len)",
      Rdist = "dCJS_ss(probSurvive, probCapture, len = 0)",
      discrete = TRUE,
      types = c('value = double(1)', 'probSurvive = double()', 'probCapture = double()', 'len = integer()'),
      pqAvail = FALSE)), verbose = F
    )

  registerDistributions(list(
    dCJS_sv = list(
      BUGSdist = "dCJS_sv(probSurvive, probCapture, len)",
      Rdist = "dCJS_sv(probSurvive, probCapture, len = 0)",
      discrete = TRUE,
      types = c('value = double(1)', 'probSurvive = double()', 'probCapture = double(1)', 'len = integer()'),
      pqAvail = FALSE)), verbose = F
    )

  registerDistributions(list(
    dCJS_vs = list(
      BUGSdist = "dCJS_vs(probSurvive, probCapture, len)",
      Rdist = "dCJS_vs(probSurvive, probCapture, len = 0)",
      discrete = TRUE,
      types = c('value = double(1)', 'probSurvive = double(1)', 'probCapture = double()', 'len = integer()'),
      pqAvail = FALSE)), verbose = F
    )

  registerDistributions(list(
    dCJS_vv = list(
      BUGSdist = "dCJS_vv(probSurvive, probCapture, len)",
      Rdist = "dCJS_vv(probSurvive, probCapture, len = 0)",
      discrete = TRUE,
      types = c('value = double(1)', 'probSurvive = double(1)', 'probCapture = double(1)', 'len = integer()'),
      pqAvail = FALSE)), verbose = F
  )

    registerDistributions(list(
      dOcc_s = list(
        BUGSdist = "dOcc_s(probOcc, probDetect, len)",
        Rdist = "dOcc_s(probOcc, probDetect, len)",
        discrete = TRUE,
        types = c('value = double(1)', 'probOcc = double(0)', 'probDetect = double(0)', 'len = integer(0)'),
        pqAvail = FALSE)), verbose = F
      )

    registerDistributions(list(
      dOcc_v = list(
        BUGSdist = "dOcc_v(probOcc, probDetect, len)",
        Rdist = c("dOcc_v(probOcc, probDetect, len)"),
        discrete = TRUE,
        types = c('value = double(1)', 'probOcc = double(0)', 'probDetect = double(1)', 'len = integer(0)'),
        pqAvail = FALSE)), verbose = F
      )


  registerDistributions(list(
      dDynOcc_vvm = list(
          BUGSdist = "dDynOcc_vvm(init, probPersist, probColonize, p, start, end)",
          Rdist = "dDynOcc_vvm(init, probPersist, probColonize, p, start, end)",
          types = c('value = double(2)',
                    'init = double()',
                    'probPersist = double(1)',
                    'probColonize = double(1)',
                    'p = double(2)',
                    'start = double(1)',
                    'end = double(1)'),
          mixedSizes = TRUE)), verbose = F
      )
  registerDistributions(list(
      dDynOcc_vsm = list(
          BUGSdist = "dDynOcc_vsm(init, probPersist, probColonize, p, start, end)",
          Rdist = "dDynOcc_vsm(init, probPersist, probColonize, p, start, end)",
          types = c('value = double(2)',
                    'init = double()',
                    'probPersist = double(1)',
                    'probColonize = double()',
                    'p = double(2)',
                    'start = double(1)',
                    'end = double(1)'),
          mixedSizes = TRUE)), verbose = F
      )

  registerDistributions(list(
      dDynOcc_svm = list(
          BUGSdist = "dDynOcc_svm(init, probPersist, probColonize, p, start, end)",
          Rdist = "dDynOcc_svm(init, probPersist, probColonize, p, start, end)",
          types = c('value = double(2)',
                    'init = double()',
                    'probPersist = double()',
                    'probColonize = double(1)',
                    'p = double(2)',
                    'start = double(1)',
                    'end = double(1)'),
          mixedSizes = TRUE)), verbose = F
      )

  registerDistributions(list(
      dDynOcc_ssm = list(
          BUGSdist = "dDynOcc_ssm(init, probPersist, probColonize, p, start, end)",
          Rdist = "dDynOcc_ssm(init, probPersist, probColonize, p, start, end)",
          types = c('value = double(2)',
                    'init = double()',
                    'probPersist = double()',
                    'probColonize = double()',
                    'p = double(2)',
                    'start = double(1)',
                    'end = double(1)'),
          mixedSizes = TRUE)), verbose = F
      )


  registerDistributions(list(
      dDynOcc_vvv = list(
          BUGSdist = "dDynOcc_vvv(init, probPersist, probColonize, p, start, end)",
          Rdist = "dDynOcc_vvv(init, probPersist, probColonize, p, start, end)",
          types = c('value = double(2)',
                    'init = double()',
                    'probPersist = double(1)',
                    'probColonize = double(1)',
                    'p = double(1)',
                    'start = double(1)',
                    'end = double(1)'),
          mixedSizes = TRUE)), verbose = F
      )
  registerDistributions(list(
      dDynOcc_vsv = list(
          BUGSdist = "dDynOcc_vsv(init, probPersist, probColonize, p, start, end)",
          Rdist = "dDynOcc_vsv(init, probPersist, probColonize, p, start, end)",
          types = c('value = double(2)',
                    'init = double()',
                    'probPersist = double(1)',
                    'probColonize = double()',
                    'p = double(1)',
                    'start = double(1)',
                    'end = double(1)'),
          mixedSizes = TRUE)), verbose = F
      )

  registerDistributions(list(
      dDynOcc_svv = list(
          BUGSdist = "dDynOcc_svv(init, probPersist, probColonize, p, start, end)",
          Rdist = "dDynOcc_svv(init, probPersist, probColonize, p, start, end)",
          types = c('value = double(2)',
                    'init = double()',
                    'probPersist = double()',
                    'probColonize = double(1)',
                    'p = double(1)',
                    'start = double(1)',
                    'end = double(1)'),
          mixedSizes = TRUE)), verbose = F
      )

  registerDistributions(list(
      dDynOcc_ssv = list(
          BUGSdist = "dDynOcc_ssv(init, probPersist, probColonize, p, start, end)",
          Rdist = "dDynOcc_ssv(init, probPersist, probColonize, p, start, end)",
          types = c('value = double(2)',
                    'init = double()',
                    'probPersist = double()',
                    'probColonize = double()',
                    'p = double(1)',
                    'start = double(1)',
                    'end = double(1)'),
          mixedSizes = TRUE)), verbose = F
      )


  registerDistributions(list(
      dDynOcc_vvs = list(
          BUGSdist = "dDynOcc_vvs(init, probPersist, probColonize, p, start, end)",
          Rdist = "dDynOcc_vvs(init, probPersist, probColonize, p, start, end)",
          types = c('value = double(2)',
                    'init = double()',
                    'probPersist = double(1)',
                    'probColonize = double(1)',
                    'p = double()',
                    'start = double(1)',
                    'end = double(1)'),
          mixedSizes = TRUE)), verbose = F
      )
  registerDistributions(list(
      dDynOcc_vss = list(
          BUGSdist = "dDynOcc_vss(init, probPersist, probColonize, p, start, end)",
          Rdist = "dDynOcc_vss(init, probPersist, probColonize, p, start, end)",
          types = c('value = double(2)',
                    'init = double()',
                    'probPersist = double(1)',
                    'probColonize = double()',
                    'p = double()',
                    'start = double(1)',
                    'end = double(1)'),
          mixedSizes = TRUE)), verbose = F
      )

  registerDistributions(list(
      dDynOcc_svs = list(
          BUGSdist = "dDynOcc_svs(init, probPersist, probColonize, p, start, end)",
          Rdist = "dDynOcc_svs(init, probPersist, probColonize, p, start, end)",
          types = c('value = double(2)',
                    'init = double()',
                    'probPersist = double()',
                    'probColonize = double(1)',
                    'p = double()',
                    'start = double(1)',
                    'end = double(1)'),
          mixedSizes = TRUE)), verbose = F
      )

  registerDistributions(list(
      dDynOcc_sss = list(
          BUGSdist = "dDynOcc_sss(init, probPersist, probColonize, p, start, end)",
          Rdist = "dDynOcc_sss(init, probPersist, probColonize, p, start, end)",
          types = c('value = double(2)',
                    'init = double()',
                    'probPersist = double()',
                    'probColonize = double()',
                    'p = double()',
                    'start = double(1)',
                    'end = double(1)'),
          mixedSizes = TRUE)), verbose = F
      )
  registerDistributions(list(
    dHMM = list(
      BUGSdist = "dHMM(init, probObs, probTrans, len, checkRowSums)",
      Rdist = "dHMM(init, probObs, probTrans, len = 0, checkRowSums = 1)",
      discrete = TRUE,
      types = c('value = double(1)',
                'init = double(1)',
                'probObs = double(2)',
                'probTrans = double(2)',
                'len = integer(0)',
                'checkRowSums = integer(0)'),
      mixedSizes = TRUE,
      pqAvail = FALSE)), verbose = F
    )

  registerDistributions(list(
    dHMMo = list(
      BUGSdist = "dHMMo(init, probObs, probTrans, len, checkRowSums)",
      Rdist = "dHMMo(init, probObs, probTrans, len = 0, checkRowSums = 1)",
      discrete = TRUE,
      types = c('value = double(1)',
                'init = double(1)',
                'probObs = double(3)',
                'probTrans = double(2)',
                'len = integer(0)',
                'checkRowSums = integer(0)'),
      mixedSizes = TRUE,
      pqAvail = FALSE)), verbose = F
    )
  registerDistributions(list(
    dDHMM = list(
      BUGSdist = "dDHMM(init, probObs, probTrans, len, checkRowSums)",
      Rdist = "dDHMM(init, probObs, probTrans, len, checkRowSums)",
      discrete = TRUE,
      types = c('value = double(1)',
                'init = double(1)',
                'probObs = double(2)',
                'probTrans = double(3)',
                'len = integer()',
                'checkRowSums = integer(0)'),
      mixedSizes = TRUE,
      pqAvail = FALSE)), verbose = F
    )

  registerDistributions(list(
    dDHMMo = list(
      BUGSdist = "dDHMMo(init, probObs, probTrans, len, checkRowSums)",
      Rdist = "dDHMMo(init, probObs, probTrans, len, checkRowSums)",
      discrete = TRUE,
      types = c('value = double(1)',
                'init = double(1)',
                'probObs = double(3)',
                'probTrans = double(3)',
                'len = integer()',
                'checkRowSums = integer(0)'),
      mixedSizes = TRUE,
      pqAvail = FALSE)), verbose = F
    )

  registerDistributions(list(
    dNmixture_v = list(
      BUGSdist = "dNmixture_v(lambda, prob, Nmin, Nmax, len)",
      Rdist = "dNmixture_v(lambda, prob, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'prob = double(1)',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
                ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )

  registerDistributions(list(
    dNmixtureAD_v = list(
      BUGSdist = "dNmixtureAD_v(lambda, prob, Nmin, Nmax, len)",
      Rdist = "dNmixtureAD_v(lambda, prob, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'prob = double(1)',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
                ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )

  registerDistributions(list(
    dNmixture_s = list(
      BUGSdist = "dNmixture_s(lambda, prob, Nmin, Nmax, len)",
      Rdist = "dNmixture_s(lambda, prob, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'prob = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
                ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )

  registerDistributions(list(
    dNmixtureAD_s = list(
      BUGSdist = "dNmixtureAD_s(lambda, prob, Nmin, Nmax, len)",
      Rdist = "dNmixtureAD_s(lambda, prob, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'prob = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
                ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )


  registerDistributions(list(
    dNmixture_BNB_v = list(
      BUGSdist = "dNmixture_BNB_v(lambda, theta, prob, Nmin, Nmax, len)",
      Rdist = "dNmixture_BNB_v(lambda, theta, prob, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'theta = double()',
                'prob = double(1)',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )
  registerDistributions(list(
    dNmixtureAD_BNB_v = list(
      BUGSdist = "dNmixtureAD_BNB_v(lambda, theta, prob, Nmin, Nmax, len)",
      Rdist = "dNmixtureAD_BNB_v(lambda, theta, prob, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'theta = double()',
                'prob = double(1)',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )
  registerDistributions(list(
    dNmixture_BNB_s = list(
      BUGSdist = "dNmixture_BNB_s(lambda, theta, prob, Nmin, Nmax, len)",
      Rdist = "dNmixture_BNB_s(lambda, theta, prob, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'theta = double()',
                'prob = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )
  registerDistributions(list(
    dNmixtureAD_BNB_s = list(
      BUGSdist = "dNmixtureAD_BNB_s(lambda, theta, prob, Nmin, Nmax, len)",
      Rdist = "dNmixtureAD_BNB_s(lambda, theta, prob, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'theta = double()',
                'prob = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )
  registerDistributions(list(
    dNmixture_BNB_oneObs = list(
      BUGSdist = "dNmixture_BNB_oneObs(lambda, theta, prob, Nmin, Nmax)",
      Rdist = "dNmixture_BNB_oneObs(lambda, theta, prob, Nmin, Nmax)",
      discrete = TRUE,
      types = c('value = double()',
                'lambda = double()',
                'theta = double()',
                'prob = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )
  registerDistributions(list(
    dNmixtureAD_BNB_oneObs = list(
      BUGSdist = "dNmixtureAD_BNB_oneObs(lambda, theta, prob, Nmin, Nmax)",
      Rdist = "dNmixtureAD_BNB_oneObs(lambda, theta, prob, Nmin, Nmax)",
      discrete = TRUE,
      types = c('value = double()',
                'lambda = double()',
                'theta = double()',
                'prob = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )

  registerDistributions(list(
    dNmixture_BBP_v = list(
      BUGSdist = "dNmixture_BBP_v(lambda, prob, s, Nmin, Nmax, len)",
      Rdist = "dNmixture_BBP_v(lambda, prob, s, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'prob = double(1)',
                's = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )
  registerDistributions(list(
    dNmixtureAD_BBP_v = list(
      BUGSdist = "dNmixtureAD_BBP_v(lambda, prob, s, Nmin, Nmax, len)",
      Rdist = "dNmixtureAD_BBP_v(lambda, prob, s, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'prob = double(1)',
                's = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )

  registerDistributions(list(
    dNmixture_BBP_s = list(
      BUGSdist = "dNmixture_BBP_s(lambda, prob, s, Nmin, Nmax, len)",
      Rdist = "dNmixture_BBP_s(lambda, prob, s, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'prob = double()',
                's = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )
  registerDistributions(list(
    dNmixtureAD_BBP_s = list(
      BUGSdist = "dNmixtureAD_BBP_s(lambda, prob, s, Nmin, Nmax, len)",
      Rdist = "dNmixtureAD_BBP_s(lambda, prob, s, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'prob = double()',
                's = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )

  registerDistributions(list(
    dNmixture_BBP_oneObs = list(
      BUGSdist = "dNmixture_BBP_oneObs(lambda, prob, s, Nmin, Nmax)",
      Rdist = "dNmixture_BBP_oneObs(lambda, prob, s, Nmin, Nmax)",
      discrete = TRUE,
      types = c('value = double()',
                'lambda = double()',
                'prob = double()',
                's = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )
  registerDistributions(list(
    dNmixtureAD_BBP_oneObs = list(
      BUGSdist = "dNmixtureAD_BBP_oneObs(lambda, prob, s, Nmin, Nmax)",
      Rdist = "dNmixtureAD_BBP_oneObs(lambda, prob, s, Nmin, Nmax)",
      discrete = TRUE,
      types = c('value = double()',
                'lambda = double()',
                'prob = double()',
                's = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )




  registerDistributions(list(
    dNmixture_BBNB_v = list(
      BUGSdist = "dNmixture_BBNB_v(lambda, theta, prob, s, Nmin, Nmax, len)",
      Rdist = "dNmixture_BBNB_v(lambda, theta, prob, s, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'prob = double(1)',
                's = double()',
                'theta = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )

  registerDistributions(list(
    dNmixtureAD_BBNB_v = list(
      BUGSdist = "dNmixtureAD_BBNB_v(lambda, theta, prob, s, Nmin, Nmax, len)",
      Rdist = "dNmixtureAD_BBNB_v(lambda, theta, prob, s, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'prob = double(1)',
                's = double()',
                'theta = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )

  registerDistributions(list(
    dNmixture_BBNB_s = list(
      BUGSdist = "dNmixture_BBNB_s(lambda, theta, prob, s, Nmin, Nmax, len)",
      Rdist = "dNmixture_BBNB_s(lambda, theta, prob, s, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'prob = double()',
                's = double()',
                'theta = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )
  registerDistributions(list(
    dNmixtureAD_BBNB_s = list(
      BUGSdist = "dNmixtureAD_BBNB_s(lambda, theta, prob, s, Nmin, Nmax, len)",
      Rdist = "dNmixtureAD_BBNB_s(lambda, theta, prob, s, Nmin, Nmax, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'lambda = double()',
                'prob = double()',
                's = double()',
                'theta = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )

    registerDistributions(list(
    dNmixture_BBNB_oneObs = list(
      BUGSdist = "dNmixture_BBNB_oneObs(lambda, theta, prob, s, Nmin, Nmax)",
      Rdist = "dNmixture_BBNB_oneObs(lambda, theta, prob, s, Nmin, Nmax)",
      discrete = TRUE,
      types = c('value = double()',
                'lambda = double()',
                'theta = double()',
                'prob = double()',
                's = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )
    registerDistributions(list(
    dNmixtureAD_BBNB_oneObs = list(
      BUGSdist = "dNmixtureAD_BBNB_oneObs(lambda, theta, prob, s, Nmin, Nmax)",
      Rdist = "dNmixtureAD_BBNB_oneObs(lambda, theta, prob, s, Nmin, Nmax)",
      discrete = TRUE,
      types = c('value = double()',
                'lambda = double()',
                'theta = double()',
                'prob = double()',
                's = double()',
                'Nmin = double(0, default = -1)',
                'Nmax = double(0, default = -1)'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )

  registerDistributions(list(
    dBetaBinom_v = list(
      BUGSdist = "dBetaBinom_v(N, shape1, shape2, len)",
      Rdist = "dBetaBinom_v(N, shape1, shape2, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'N = double()',
                'shape1 = double(1)',
                'shape2 = double(1)',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )

  registerDistributions(list(
    dBetaBinom_s = list(
      BUGSdist = "dBetaBinom_s(N, shape1, shape2, len)",
      Rdist = "dBetaBinom_s(N, shape1, shape2, len)",
      discrete = TRUE,
      types = c('value = double(1)',
                'N = double()',
                'shape1 = double()',
                'shape2 = double()',
                'len = double()'
      ),
      mixedSizes = FALSE,
      pqAvail = FALSE
    )), verbose = F
  )


})}
