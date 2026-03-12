#' Simulated Dataset Generation for High-Dimensional Mediation Analysis
#'
#' Generates synthetic datasets mimicking high-dimensional mediation structures,
#' optionally incorporating transferable source data under varying covariate correlation
#' and heterogeneity levels. This function supports heterogeneous settings for data
#' generation.
#'
#' @param n Integer. Number of observations (sample size). Default is 100.
#' @param p_x Integer. Number of covariates (confounders). Default is 5.
#' @param rho Numeric. Correlation coefficient (0–1) controlling correlation between mediators. Default is 0 (no correlation).
#' @param p_m Integer. Number of mediators. Default is 100.
#' @param h Integer. Degree of heterogeneity (for source data). Default is 0.
#' @param source Logical. If TRUE, generate source (external) dataset. Default is FALSE.
#' @param transferable Logical. If TRUE, generates a transferable source dataset sharing mediator-outcome structure with the target. Default is TRUE.
#' @param seed Integer. Optional random seed for reproducibility. Default is NULL.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item{\code{data}}: A data.frame of dimension \code{n × (2 + p_m + p_x)} containing outcome (Y), treatment (D), mediators (M1–Mp_m), and covariates (X1–Xp_x).
#'   \item{\code{coef}}: A named list of true model coefficients (\code{alpha1, alpha2, beta1, beta2, beta4}, etc.).
#' }
#'
#' @details
#' This function generates data according to a structural equation model (SEM):
#' \itemize{
#'   \item Treatment D depends linearly on covariates X.
#'   \item Mediators M depend on D and X, with residual correlation controlled by \code{rho}.
#'   \item Outcome Y depends on D, M, and X.
#'   \item If \code{source = TRUE}, then a source dataset is simulated with potentially transferable mechanisms.
#' }
#'
#' @examples
#' source_data <- gen_simData_hetero(
#'   n = 100, p_x = 5, rho = 0, p_m = 100, h = 0,
#'   source = TRUE, transferable = TRUE, seed = 123
#' )
#' source_data <- gen_simData_hetero(
#'   n = 100, p_x = 5, rho = 0, p_m = 100, h = 0,
#'   source = FALSE, transferable = TRUE, seed = 123
#' )
#' @export
gen_simData_hetero<-function(n=100,
                             p_x=5,
                             rho=0,
                             p_m=100,
                             h=0,
                             source=FALSE,
                             transferable=TRUE,
                             seed=NULL){

  if(!is.null(seed)){
    set.seed(seed)
  }
  gamma<-rep(0.3,p_x)
  alpha0<-0
  alpha1<-c(c(0.6,0.55,0.65,0.5,0.35,0.3),0,0,rep(0,p_m-8))
  alpha2<-sapply(1:p_m,function(x){rep(0.3,p_x)})

  beta0<-0
  beta1<-1
  beta2<-c(c(0.4,0.45,0.45,0.5),0,0,c(0.35,0.25),rep(0,p_m-8))
  beta3<-rep(0,p_m)
  beta4<-rep(0.3,p_x)

  if(source){
    if(transferable){
      alpha1<-c(c(0.6,0.55,0.65,0.5,0.35,0.3),0,0,sample(c(rep(0.3,h),rep(0,p_m-8-h))))
      beta2<-c(c(0.4,0.45,0.45,0.5),0,0,c(0.35,0.25),sample(c(rep(0.3,h),rep(0,p_m-8-h))))
    }else{
      alpha1<-c(rep(0.5,8),sample(c(rep(0.5,2),rep(0,p_m-8-2))))
      beta2<-c(rep(0.5,8),sample(c(rep(0.5,2),rep(0,p_m-8-2))))
    }
  }

  ### confounder: X
  X.Sigma=diag(rep(1,p_x))
  for(i in 1:p_x){for(j in 1:p_x){X.Sigma[i,j]=0.5^abs(i-j)}}

  if(source){
    ### covariate shift
    A_k <- matrix(
      data = sample(c(0, 0.3), size = p_x * p_x, replace = TRUE, prob = c(0.7, 0.3)),
      nrow = p_x,
      ncol = p_x
    )

    X.Sigma<-X.Sigma+t(A_k)%*%A_k
  }

  X<-MASS::mvrnorm(n, c(rep(0,p_x)), X.Sigma)

  ### treatment: D
  D <- X%*%gamma + rnorm(n,0,1)

  ### mediator
  M.Sigma=diag(rep(1,p_m))
  for(i in 1:p_m){for(j in 1:p_m){M.Sigma[i,j]=rho^abs(i-j)}}
  M.se<-MASS::mvrnorm(n,rep(0,p_m), M.Sigma)
  M<-matrix(data=NA,nrow=n,ncol=p_m)
  for(i in 1:p_m){
    M[,i] <- alpha1[i]*D+X%*%alpha2[,i]+M.se[,i]
  }

  ### outcome
  Y<-beta0 + beta1*D + M%*%beta2 + X%*%beta4 + rnorm(n,0,1)

  data<-as.data.frame(cbind(Y=Y,D=D,M=M,X=X))
  colnames(data)<-c('Y','D',paste0('M',1:p_m),paste0('X',1:p_x))

  ### true effect
  coef<-list(alpha0=unname(alpha0),
             alpha1=unname(alpha1),
             alpha2=unname(alpha2),
             beta0=unname(beta0),
             beta1=unname(beta1),
             beta2=unname(beta2),
             beta3=unname(beta3),
             beta4=unname(beta4))

  return(list(data=data,coef=coef))
}
