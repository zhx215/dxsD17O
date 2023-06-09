# developed by ZX and MJW, 2021
# stand-alone function for calculating raindrop re-evaporation

source('lcl.R')
reevap <- function(p0,temp0,rh0,D,d18Ov,d17Ov,d2Hv,d2H_offset=0,d17O_offset=0) { # inputs: surface pressure, surface temperature, surface RH, initial raindrop diameter, vapor d18O, vapor d2H 
  a18eff <- function(temp,mode="k",l=0.004) {
    if (temp >= 273.15) {
      eq <- exp(1137/temp^2-0.4156/temp-0.00207)
      value <- eq
    }
    else if (temp <= 250.15) {
      eq <- exp(11.839/temp-0.028224)
      Si <- 1-l*(temp-273.15)
      diff16_18 <- 1.0285
      k <- Si/(eq*diff16_18*(Si-1)+1)
      if (mode=="no k") {k <- 1}
      value <- k*eq
    }
    else {
      eq_273 <- exp(1137/273.15^2-0.4156/273.15-0.00207)
      eq_274 <- exp(1137/274.15^2-0.4156/274.15-0.00207)
      eq_250 <- exp(11.839/250.15-0.028224)
      eq_249 <- exp(11.839/249.15-0.028224)
      Si_250 <- 1-l*(250.15-273.15)
      Si_249 <- 1-l*(249.15-273.15)
      diff16_18 <- 1.0285
      k_250 <- Si_250/(eq_250*diff16_18*(Si_250-1)+1)
      k_249 <- Si_249/(eq_249*diff16_18*(Si_249-1)+1)
      if (mode=="no k") {k_250 <- k_249 <- 1}
      x1 <- 250.15
      x2 <- 273.15
      k1 <- eq_250*k_250-eq_249*k_249
      k2 <- eq_274-eq_273
      y1 <- eq_250*k_250
      y2 <- eq_273
      answer <- solve(matrix(c(x1^3,x2^3,3*x1^2,3*x2^2,x1^2,x2^2,2*x1,2*x2,x1,x2,1,1,1,1,0,0),ncol=4),c(y1,y2,k1,k2))
      value <- answer[1]*temp^3+answer[2]*temp^2+answer[3]*temp+answer[4]
    }
    return(value)
  }
  
  a17eff <- function(temp,mode="k",l=0.004) {
    if (temp >= 273.15) {
      eq <- exp(1137/temp^2-0.4156/temp-0.00207)^0.529
      value <- eq
    }
    else if (temp <= 250.15) {
      eq <- exp(11.839/temp-0.028224)^0.529
      Si <- 1-l*(temp-273.15)
      diff16_17 <- 1.0285^0.518
      k <- Si/(eq*diff16_17*(Si-1)+1)
      if (mode=="no k") {k <- 1}
      value <- k*eq
    }
    else {
      eq_273 <- exp(1137/273.15^2-0.4156/273.15-0.00207)^0.529
      eq_274 <- exp(1137/274.15^2-0.4156/274.15-0.00207)^0.529
      eq_250 <- exp(11.839/250.15-0.028224)^0.529
      eq_249 <- exp(11.839/249.15-0.028224)^0.529
      Si_250 <- 1-l*(250.15-273.15)
      Si_249 <- 1-l*(249.15-273.15)
      diff16_17 <- 1.0285^0.518
      k_250 <- Si_250/(eq_250*diff16_17*(Si_250-1)+1)
      k_249 <- Si_249/(eq_249*diff16_17*(Si_249-1)+1)
      if (mode=="no k") {k_250 <- k_249 <- 1}
      x1 <- 250.15
      x2 <- 273.15
      k1 <- eq_250*k_250-eq_249*k_249
      k2 <- eq_274-eq_273
      y1 <- eq_250*k_250
      y2 <- eq_273
      answer <- solve(matrix(c(x1^3,x2^3,3*x1^2,3*x2^2,x1^2,x2^2,2*x1,2*x2,x1,x2,1,1,1,1,0,0),ncol=4),c(y1,y2,k1,k2))
      value <- answer[1]*temp^3+answer[2]*temp^2+answer[3]*temp+answer[4]
    }
    return(value)
  }
  
  a2eff <- function(temp,mode="k",l=0.004) {
    if (temp >= 273.15) {
      eq <- exp(24844/temp^2-76.248/temp+0.05261)
      value <- eq
    }
    else if (temp <= 250.15) {
      eq <- exp(16289/temp^2-0.0945)
      Si <- 1-l*(temp-273.15)
      diff1_2 <- 1.0251
      k <- Si/(eq*diff1_2*(Si-1)+1)
      if (mode=="no k") {k <- 1}
      value <- k*eq
    }
    else {
      eq_273 <- exp(24844/273.15^2-76.248/273.15+0.05261)
      eq_274 <- exp(24844/274.15^2-76.248/274.15+0.05261)
      eq_250 <- exp(16289/250.15^2-0.0945)
      eq_249 <- exp(16289/249.15^2-0.0945)
      Si_250 <- 1-l*(250.15-273.15)
      Si_249 <- 1-l*(249.15-273.15)
      diff1_2 <- 1.0251
      k_250 <- Si_250/(eq_250*diff1_2*(Si_250-1)+1)
      k_249 <- Si_249/(eq_249*diff1_2*(Si_249-1)+1)
      if (mode=="no k") {k_250 <- k_249 <- 1}
      x1 <- 250.15
      x2 <- 273.15
      k1 <- eq_250*k_250-eq_249*k_249
      k2 <- eq_274-eq_273
      y1 <- eq_250*k_250
      y2 <- eq_273
      answer <- solve(matrix(c(x1^3,x2^3,3*x1^2,3*x2^2,x1^2,x2^2,2*x1,2*x2,x1,x2,1,1,1,1,0,0),ncol=4),c(y1,y2,k1,k2))
      value <- answer[1]*temp^3+answer[2]*temp^2+answer[3]*temp+answer[4]
    }
    return(value)
  }
  
  # if RH = 1, no raindrop re-evaporation as LCL height = 0 m
  if (rh0 == 1){
    return(c(0,1,1,1,D))
  }
  
  # saturation vapor pressure over liquid (Romps, 2017)
  get_es <- function(temp){ # K
    Ptrip <- 611.65 # Pa
    Ttrip <- 273.16 # K
    E0v <- 2.374*10^6 # J/kg
    cvl <- 4119 # J/kg/K, specific heat capacity of liquid water
    cvv <- 1418 # J/kg/K, specific heat capacity of water vapor at constant volume
    Rv <- 461 # J/kg/K, specific gas constant of water vapor
    cpv <- cvv + Rv # J/kg/K, specific heat capacity of water vapor at constant pressure
    Ptrip*(temp/Ttrip)^((cpv-cvl)/Rv)*exp(((E0v-(cvv-cvl)*Ttrip)/Rv)*(1/Ttrip-1/temp)) # Pa
  }
  
  # pressure at a temperature of certain level, given initial P/T/RH
  get_p <- function(p0,temp0,rh0,t_z){
    Rv <- 461 # J/kg/K, specific gas constant of water vapor
    Ra <- 287.04 # J/kg/K, specific gas constant of dry air
    pv <- rh0*get_es(temp0) # Pa, partial pressure of water vapor
    qv <- Ra*pv/(Rv*p0+(Ra-Rv)*pv) # mass fraction of water vapor
    Rm <- (1-qv)*Ra+qv*Rv # J/kg/K, air parcel's specific gas constant
    cvv <- 1418 # J/kg/K, specific heat capacity of water vapor at constant volume
    cpv <- cvv + Rv # J/kg/K, specific heat capacity of water vapor at constant pressure
    cva <- 719 # J/kg/K, specific heat capacity of dry air at constant volume
    cpa <- cva + Ra # J/kg/K, heat capacity of dry air at constant pressure
    cpm <- (1-qv)*cpa + qv*cpv # J/kg/K, air parcel's specific heat capacity at constant pressure
    p0*(t_z/temp0)^(cpm/Rm)
  }
  
  # get pressure at a temperature of certain level, given initial P/T/RH
  get_rh <- function(p0,temp0,rh0,t_z){
    es_c <- get_es(lcl_data[2])
    p_c <- get_p(p0,temp0,rh0,lcl_data[2])
    ws_c <- 0.622*es_c/(p_c-es_c)
    es <- get_es(t_z)
    p <- get_p(p0,temp0,rh0,t_z)
    ws <-0.622*es/(p-es)
    es0 <- get_es(temp0)
    ws0 <- 0.622*es0/(p0-es0)
    w <- ws0*rh0
    rh <- w/ws
    rh_c <- w/ws_c
    (rh-rh0)/(rh_c-rh0)*(1-rh0)+rh0 # normalize final RH to 0-1 scale
  }
  
  # get effective RH for raindrops
  get_rheff <- function(temp,temp_d,rh){
    rheff <- get_es(temp)*rh/get_es(temp_d)
    if (rheff>1) {
      return(1) # set rheff as 1 if > 1
    } else {
      return(rheff)
    }
  }
  
  # air density, Salamalikis et al., 2016, Plathner and Woloszyn, 2002
  rho <- function(p,temp,rh,es){
    Rd <- 287.04 # J/kg/K, specific gas constant of dry air
    1/(Rd*temp)*(p-0.378*rh*es) # kg/m3
  }
  
  # terminal velocity, Foote and Du Toit, 1969
  tv <- function(d,rho,rho_ref){
    9.43*(1-exp(-(d/1.77)^1.147))*(rho_ref/rho)^0.4 # m/s 
  }
  
  # dynamic viscosity of air, Sutherland, 1893
  visc <- function(temp){
    1.458*10^(-6)*temp^1.5/(temp+110.4) # kg/m/s
  }
  
  # diffusivity of water molecules, Hall and Pruppacher, 1976
  diff <- function(p,temp){
    p0 <- 101325
    0.211*10^-4*(p0/p)*(temp/273.15)^1.94 # m2/s
  }
  
  # water density, based on drop temperature, slightly offset from 1000 kg/m3 or 1 g/cm3, Salamalikis et al., 2016
  rhow <- function(temp_d){
    t <- temp_d-273.15
    1000*(1-(t+288.9414)*(t-3.9863)^2/508929.2/(t+68.12963)) # kg/m3
  }
  
  # latent heat of evaporation, based on drop temperature, Rogers and Yau, 1989
  Le <- function(temp_d){
    t <- temp_d-273.15
    2500800-2360*t+1.6*t^2-0.06*t^3 # J/kg
  }
  
  # thermal conductivity of air, ka is close to k, Pruppacher and Klett, 2010
  ka <- function(temp){
    418.4*(5.69+0.017*(temp-273.15))*10^-5 # J/m/s/K, the original form is in unit cal/cm/s/K (1 cal = 4.184 J)
  }
  
  # initalize the model
  lcl_data <- lcl(p0,temp0,rhl=rh0) # get lcl data
  gamma <- (temp0-lcl_data[2])/lcl_data[1] # K/m, lapse rate
  z_p <- lcl_data[1] # height profile
  D_p <- D # drop diameter profile
  m_p <- 4/3*pi*(D_p/2)^3*10^-3 # g, mass profile
  t_p <- td_p <- lcl_data[2] # temperature and drop temperature profile
  p_p <- get_p(p0,temp0,rh0,t_p)# pressure profile
  rh_p <- get_rh(p0,temp0,rh0,t_p) # rh profile
  rheff_p <- get_rheff(t_p,td_p,rh_p) # effective rh profile
  rho_p <- rho(p_p,t_p,rh_p,get_es(t_p)) # air density profile
  rho_ref <- rho(p0,temp0,rh0,get_es(temp0)) # reference air density
  tv_p <- tv(D_p,rho_p,rho_ref) # terminal velocity profile
  visc_p <- visc(t_p) # viscosity profile
  diff_p <- diff(p_p,t_p) # diffusivity profile
  rhow_p <- rhow(td_p) # water density profile
  Le_p <- Le(td_p) # latent heat profile
  ka_p <- ka(t_p) # conductivity profile
  Sc_p <- visc_p/(rho_p*diff_p) # Schmidt number profile
  cp <- 1004 # J/kg/K specific heat of dry air
  Pr_p <- visc_p*cp/ka_p # Prandtl number profile
  Re_p <- tv_p*(D_p/1000)*rho_p/visc_p # Reynolds number
  fv_p <- fh_p <- vector() # mass and heat ventilation coefficient
  fv_p18 <- vector() # mass ventilation coefficient for 18O
  fv_p17 <- vector()
  fv_p2 <- vector() # mass ventilation coefficient for 2H
  Rw <- 461 # J/kg/K, specific gas constant of water vapor
  cw <- 4119 # J/kg/K specific heat capacity of liquid water, Romps 2017
  diff16_18 <- 1.0285 # (D/D') for 16O/18O
  diff16_17 <- 1.0285^0.518 # (D/D') for 16O/17O
  diff1_2 <- 1.0251 # (D/D') for 1H/2H
  rstd18 <- 2005.2*10^-6 # VSMOW
  rstd17 <- 379.9*10^-6
  rstd2 <- 155.76*10^-6 # VSMOW
  rv18 <- (d18Ov/1000+1)*rstd18 # get r value for 18O/16O in vapor, could be revised to study equilibration
  rv17 <- (d17Ov/1000+1)*rstd17
  rv2 <- (d2Hv/1000+1)*rstd2 # get r value for 2H/1H in vapor, could be revised to study equilibration
  rp18_p <- rv18*a18eff(td_p) # get r value for 18O/16O in raindrop, initial value
  rp17_p <- rv17*a17eff(td_p)
  rp2_p <- rv2*a2eff(td_p) # get r value for 2H/1H in raindrop, initial value
  m18_p <- m_p*rp18_p # get 18O mass
  m17_p <- m_p*rp17_p
  m2_p <- m_p*rp2_p # get 2H mass
  i <- 1
  
  rv2 <- ((d2Hv+d2H_offset)/1000+1)*rstd2
  rv17 <- ((d17Ov+d17O_offset)/1000+1)*rstd17
  
  # iterate the model
  repeat{
    # mass ventilation coefficient
    if ((Sc_p[i]^(1/3)*Re_p[i]^0.5)>=1.4){
      fv_p[i] <- 0.78+0.308*Sc_p[i]^(1/3)*Re_p[i]^0.5
      fv_p18[i] <- 0.78+0.308*(Sc_p[i]*(diff16_18))^(1/3)*Re_p[i]^0.5
      fv_p17[i] <- 0.78+0.308*(Sc_p[i]*(diff16_17))^(1/3)*Re_p[i]^0.5
      fv_p2[i] <- 0.78+0.308*(Sc_p[i]*(diff1_2))^(1/3)*Re_p[i]^0.5
    } else{
      fv_p[i] <- 1+0.108*(Sc_p[i]^(1/3)*Re_p[i]^0.5)^2
      fv_p18[i] <- 1+0.108*((Sc_p[i]*(diff16_18))^(1/3)*Re_p[i]^0.5)^2
      fv_p17[i] <- 1+0.108*((Sc_p[i]*(diff16_17))^(1/3)*Re_p[i]^0.5)^2
      fv_p2[i] <- 1+0.108*((Sc_p[i]*(diff1_2))^(1/3)*Re_p[i]^0.5)^2
    }
    
    # heat ventilation coefficient
    if ((Pr_p[i]^(1/3)*Re_p[i]^0.5)>=1.4){
      fh_p[i] <- 0.78+0.308*Pr_p[i]^(1/3)*Re_p[i]^0.5
    } else{
      fh_p[i] <- 1+0.108*(Pr_p[i]^(1/3)*Re_p[i]^0.5)^2
    }
    
    # mass and heat transfer equations
    if (t_p[i]<=273.15){ # if surrounding air is frozen, no evaporation
      e_rate <- e_rate_18 <- e_rate_17 <- e_rate_2 <- 0
    } else {
      e_rate <- 2*pi*(D_p[i]/1000)*fv_p[i]*diff_p[i]/Rw*(rh_p[i]*get_es(t_p[i])/t_p[i]-get_es(td_p[i])/td_p[i]) # kg/s, evaporation rate
      e_rate_18 <- 2*pi*(D_p[i]/1000)*fv_p[i]*diff_p[i]/Rw*(fv_p18[i]/fv_p[i]/diff16_18)^0.58*
        (rh_p[i]*rv18*get_es(t_p[i])/t_p[i]-rp18_p[i]/a18eff(td_p[i],"no k")*get_es(td_p[i])/td_p[i])
      e_rate_17 <- 2*pi*(D_p[i]/1000)*fv_p[i]*diff_p[i]/Rw*(fv_p17[i]/fv_p[i]/diff16_17)^0.58*
        (rh_p[i]*rv17*get_es(t_p[i])/t_p[i]-rp17_p[i]/a17eff(td_p[i],"no k")*get_es(td_p[i])/td_p[i])
      e_rate_2 <- 2*pi*(D_p[i]/1000)*fv_p[i]*diff_p[i]/Rw*(fv_p2[i]/fv_p[i]/diff1_2)^0.58*
        (rh_p[i]*rv2*get_es(t_p[i])/t_p[i]-rp2_p[i]/a2eff(td_p[i],"no k")*get_es(td_p[i])/td_p[i])
    }
    td_rate <- 12/((D_p[i]/1000)^2*rhow_p[i]*cw)*
      (Le_p[i]*fv_p[i]*diff_p[i]/Rw*(rh_p[i]*get_es(t_p[i])/t_p[i]-get_es(td_p[i])/td_p[i])-fh_p[i]*ka_p[i]*(td_p[i]-t_p[i])) # K/s
    
    # mass and heat changes
    t <- 0.1 # s
    de <- (e_rate*1000)*t # g
    de18 <- (e_rate_18*1000)*t
    de17 <- (e_rate_17*1000)*t
    de2 <- (e_rate_2*1000)*t
    dt <- td_rate*t # K
    
    # return results if reaching the ground
    if ((z_p[i]-tv_p[i]*t)<0){
      print('hit ground')
      return(c(100-100*min(m_p)/max(m_p),tail(rp18_p,n=1)/rp18_p[1],tail(rp17_p,n=1)/rp17_p[1],tail(rp2_p,n=1)/rp2_p[1],tail(D_p,n=1))) # evaporation amount (%), alpha_rr for 18O, alpha_rr for 17O, alpha_rr for 2H, final raindrop diameter (mm)
    }
    if (abs(dt)>10^(-2)){
      print('all gone (numerical model drift)')
      return(c(100,1,1,1,0))
    }
    if ((m_p[i]+de)<0){
      print('all gone (ideally)')
      return(c(100,1,1,1,0))
    }
    
    # update model parameters
    z_p[i+1] <- z_p[i]-tv_p[i]*t
    m_p[i+1] <- m_p[i]+de
    D_p[i+1] <- (6*m_p[i+1]*10^3/pi)^(1/3)
    t_p[i+1] <- temp0-gamma*z_p[i+1]
    td_p[i+1] <- td_p[i]+dt
    p_p[i+1] <- get_p(p0,temp0,rh0,t_p[i+1])
    rh_p[i+1] <- get_rh(p0,temp0,rh0,t_p[i+1])
    rheff_p[i+1] <- get_rheff(t_p[i+1],td_p[i+1],rh_p[i+1])
    rho_p[i+1] <- rho(p_p[i+1],t_p[i+1],rh_p[i+1],get_es(t_p[i+1]))
    tv_p[i+1] <- tv(D_p[i+1],rho_p[i+1],rho_ref)
    visc_p[i+1] <- visc(t_p[i+1])
    diff_p[i+1] <- diff(p_p[i+1],t_p[i+1])
    rhow_p[i+1] <- rhow(td_p[i+1])
    Le_p[i+1] <- Le(td_p[i+1])
    ka_p[i+1] <- ka(t_p[i+1])
    Sc_p[i+1] <- visc_p[i+1]/(rho_p[i+1]*diff_p[i+1])
    Pr_p[i+1] <- visc_p[i+1]*cp/ka_p[i+1]
    Re_p[i+1] <- tv_p[i+1]*(D_p[i+1]/1000)*rho_p[i+1]/visc_p[i+1]
    m18_p[i+1] <- m18_p[i]+de18
    m17_p[i+1] <- m17_p[i]+de17
    m2_p[i+1] <- m2_p[i]+de2
    rp18_p[i+1] <- m18_p[i+1]/m_p[i+1]
    rp17_p[i+1] <- m17_p[i+1]/m_p[i+1]
    rp2_p[i+1] <- m2_p[i+1]/(m_p[i+1])
    i <- i+1
  }
}