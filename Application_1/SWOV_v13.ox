#include <oxstd.h>
#include <oxdraw.h>
#include <packages/ssfpack/ssfpack_ex.h>
#import <maximize>
#include <oxfloat.h>

static decl s_mStsm,
// variance matrix settings
s_iLevel = 2, // 0 == zeros, 1 == ones, 2 == one factor, 3 == diagonal 4 == full rank
s_iSlope = 0, // 0,1,2
s_iSeaso = 0, // 0,1,2
s_iIrreg = 4; // 1,2,3,4;
static decl s_mY, s_iN, s_iT, s_iNTstar, s_vPar, s_dVar, s_dLik, s_cst, s_dAIC, s_dBIC, s_cMod, s_vVarCorrir;

static decl s_igr = 0, s_sNameY, s_iFirstYear, s_iFirstMonth;

GetVarianceMatrixOnes(const cy, const vpar, const pipar)
{
	pipar[0] = 1;
	if (vpar == <>)
		return <>;
	else
		return exp(2.0 * vpar[0]) * ones(cy,cy);
}
GetVarianceMatrixDiagonal(const cy, const vpar, const pipar)
{
	pipar[0] = cy;
	if (vpar == <>)
		return <>;
	else
		return diag(exp(2.0 * vpar[:cy-1]));
}
GetVarianceMatrixChol(const cy, const crank, const vpar, const pipar)
{
	decl i, j, k, mA=unit(cy)[][:crank-1], mD=zeros(crank,crank);
	i = (cy-1)*cy; j = cy - crank; k = (j-1)*j;
	pipar[0] = crank + ((i-k) / 2);
	if (vpar==<>)
		return <>;
	else
	{
		for (i=0; i<crank; i++) mD[i][i] = exp(2.0 * vpar[i]);
		for (i=k=0; i<crank; i++)
			for (j=i+1; j<cy; j++, k++)
				mA[j][i] = vpar[crank+k];
		return mA * mD * mA';
	}
}

SetStsmModel1(const cy, const vpar, const fprint) // MODEL 1
{
	decl i,k,l;
	// vpar = variance level, variance slope, variance seasonal, variance matrix irreg
	s_mStsm = <CMP_LEVEL,     0,  0, 0;
		       CMP_SLOPE,     0,  0, 0;
	   		   CMP_SEAS_TRIG, 0, 12, 0;	// 12 for monthly data
			   CMP_IRREG,     1,  0, 0>;
	k = 0;
	if (s_iLevel != 0) s_mStsm[0][1] = exp(vpar[k++]); // level
	if (s_iSlope != 0) s_mStsm[1][1] = exp(vpar[k++]); // islope 1 == ones var matrix
	if (s_iSeaso != 0) s_mStsm[2][1] = exp(vpar[k++]); // iseaso 1 == ones var matrix

    decl mphi, momega, msigma, cst;
	GetSsfStsm(s_mStsm, &mphi, &momega, &msigma); // get state space model
	//if (fprint)	println("Univariate \n   - system matrices: ", mphi, "   - hyperparameters: ", momega, "   - prior: ", msigma);
	
	cst = columns(mphi);
	//if (fprint)	println("Univariate \n   - system matrices: ", mphi, "   - hyperparameters: ", momega, "   - prior: ", msigma);

	mphi = mphi ** unit(cy);
	msigma = msigma[:cst-1][:cst-1] ** unit(cy);
	msigma |= 0;
	momega = momega ** ones(cy,cy);
	decl m = zeros(cy,cy);
	i = 0;
	if (s_iIrreg == 1)
	{
		m = constant(exp(2.0*vpar[k]), cy, cy);
		i = 1;
	}
	else
	if (s_iIrreg == 2)
		m = GetVarianceMatrixChol(cy, 1, vpar[k:], &i);
	else
	if (s_iIrreg == 3)
		m = GetVarianceMatrixDiagonal(cy, vpar[k:], &i);
	else
	if (s_iIrreg == 4)
		m = GetVarianceMatrixChol(cy, cy, vpar[k:], &i);
	momega[cy*cst:][cy*cst:] = m;

	m = zeros(cy,cy);
	k += i;
	i = 0;
	l = 0;
	if (s_iLevel == 1)
		m = constant(exp(2.0*vpar[l++]), cy, cy);
	else
	if (s_iLevel == 2)
	{
		if (fprint)
		{
			println(" .. variance level = ", exp(2.0* vpar[l]));
			println(" .. loadings level = ", (1 | vpar[k:k+cy-2]));
		}

		m = GetVarianceMatrixChol(cy, 1, vpar[l++] | vpar[k:], &i);
		i--;
	}
	else
	if (s_iLevel == 3)
	{
		m = GetVarianceMatrixDiagonal(cy, vpar[l++] | vpar[k:], &i);
		i--;
	}
	else
	if (s_iLevel == 4)
	{
		m = GetVarianceMatrixChol(cy, cy, vpar[l++] | vpar[k:], &i);
		i--;
	}
	momega[:cy-1][:cy-1] = m;

	k += i;
	i=0;
	m = zeros(cy,cy);
	if (s_iSlope == 1)
		m = constant(exp(2.0*vpar[l++]), cy, cy);
	else
	if (s_iSlope == 2)
	{
		if (fprint)
		{
			println(" .. variance slope = ", exp(2.0* vpar[l]));
			println(" .. loadings slope = ", (1 | vpar[k:k+cy-2]));
		}
		m = GetVarianceMatrixChol(cy, 1, vpar[l++] | vpar[k:], &i);
		i--;
	}
	momega[cy:cy+cy-1][cy:cy+cy-1] = m;

	k += i;
	i=0;
	m = zeros(cy,cy);
	if (s_iSeaso == 1)
		m = constant(exp(2.0*vpar[l++]), cy, cy);
	else
	if (s_iSeaso == 2)
	{
		if (fprint)
		{
			println(" .. variance seasonal = ", exp(2.0* vpar[l]));
			println(" .. loadings seasonal = ", (1 | vpar[k:k+cy-2]));
		}
		m = GetVarianceMatrixChol(cy, 1, vpar[l++] | vpar[k:], &i);
		i--;
	}
	momega[cy+cy:cy+cy+(11*cy)-1][cy+cy:cy+cy+(11*cy)-1] = unit(11) ** m;
	
	// placing irregular in state vector
	cst = columns(mphi);
	decl csy = rows(mphi);   

	m = mphi[cst:][]; // Z_t
	mphi = diagcat(mphi[:cst-1][:cst-1],zeros(cy,cy)); // diag(T_t,0)
	mphi |= m ~ unit(cy); 
	m = momega[cst:][cst:];	// variance matrix irregular
	momega = diagcat(momega, zeros(cy,cy)); // diag(momega,zero matrix), meas eq has no noise
	msigma = diagcat(msigma[:cst-1][:cst-1], m) | 0; // diag(msigma,var_irreg)
	cst = columns(mphi);
	csy = rows(mphi);
	if (fprint)
	{
		println("s_iLevel: ", s_iLevel); println("s_iSlope: ", s_iSlope); println("s_iSeaso: ", s_iSeaso); println("s_iIrreg: ", s_iIrreg);
		println("vpar", vpar);
		println("mOmega level", momega[:cy-1][:cy-1]);
		decl corrmlevel, sdilevel;
		sdilevel = 1 ./ sqrt(diagonal(momega[:cy-1][:cy-1]));
		corrmlevel = diag(sdilevel) * momega[:cy-1][:cy-1] * diag(sdilevel);
		println("mOmega level corr mat", corrmlevel);
		println("mOmega slope", momega[cy:cy+cy-1][cy:cy+cy-1]);
		decl corrmslope, sdislope;
		sdislope = 1 ./ sqrt(diagonal(momega[cy:cy+cy-1][cy:cy+cy-1]));
		corrmslope = diag(sdislope) * momega[cy:cy+cy-1][cy:cy+cy-1] * diag(sdislope);
		println("mOmega slope corr mat", corrmslope);
		println("mOmega seaso", momega[cy+cy:cy+cy+cy-1][cy+cy:cy+cy+cy-1]);
		decl corrmseaso, sdiseaso;
		sdiseaso = 1 ./ sqrt(diagonal(momega[cy+cy:cy+cy+cy-1][cy+cy:cy+cy+cy-1]));
		corrmseaso = diag(sdiseaso) * momega[cy+cy:cy+cy+cy-1][cy+cy:cy+cy+cy-1] * diag(sdiseaso);
		println("mOmega seaso corr mat", corrmseaso);
		println("mOmega irreg", momega[cst-cy:cst-1][cst-cy:cst-1]);
		decl virregvar, mirregcorr, sdiirreg;
		virregvar = diagonal(momega[cst-cy:cst-1][cst-cy:cst-1]);
		println("diagonal omega irreg", virregvar);
		sdiirreg = 1 ./ sqrt(diagonal(momega[cst-cy:cst-1][cst-cy:cst-1]));
		mirregcorr = diag(sdiirreg) * momega[cst-cy:cst-1][cst-cy:cst-1] * diag(sdiirreg);
		println("mOmega irreg corr mat", mirregcorr);
		s_vVarCorrir = virregvar' | mirregcorr[1:3][0] | mirregcorr[2:3][1] | mirregcorr[3][2];
	}
	//if (fprint) println("Multivariate \n- system matrices: ", mphi, "- hyperparameters: ", momega, "- prior: ", msigma);

	
	return {mphi, momega, msigma};
}
InitPar1()
{
	decl i, j, k, l, cp;
	decl m;

	// nrpars for Irregular
	if (s_iIrreg == 1)
		cp = 1;
	else
	if (s_iIrreg == 2)
		m = GetVarianceMatrixChol(s_iN, 1, <>, &cp);
	else
	if (s_iIrreg == 3)
		m = GetVarianceMatrixDiagonal(s_iN, <>, &cp);
	else
	if (s_iIrreg == 4)
		m = GetVarianceMatrixChol(s_iN, s_iN, <>, &cp);
	
	cp += 3; // nrpar for variances of level, slope, seasonal
	if (s_iLevel == 0) cp--;
	if (s_iSlope == 0) cp--;
	if (s_iSeaso == 0) cp--;
	j = cp;
	i = 0;
	if (s_iLevel == 2)
	{
		m = GetVarianceMatrixChol(s_iN, 1, <>, &i);
		i--;
	}
	else
	if (s_iLevel == 3)
	{
		m = GetVarianceMatrixDiagonal(s_iN, <>, &i);
		i--;
	}
	else
	if (s_iLevel == 4)
	{
		m = GetVarianceMatrixChol(s_iN, s_iN, <>, &i);
		i--;
	}
	//
	k=0;
	if (s_iSlope == 2)
	{
		m = GetVarianceMatrixChol(s_iN, 1, <>, &k);
		k--; // because the +=4 for the variances was already added.
	}
	//
	l=0;
	if (s_iSeaso == 2)
	{
		m = GetVarianceMatrixChol(s_iN, 1, <>, &l);
		l--; // because the +=4 for the variances was already added.
	}
	
	cp += i+k+l;

	decl vp = zeros(cp,1);
	i=0;
	if (s_iLevel != 0) vp[i++] = -0.5; // log standard deviation level
	if (s_iSlope != 0) vp[i++] = -1.5; // log standard deviation slope 
	if (s_iSeaso != 0) vp[i++] = -1.0; // log standard deviation seasonal
	if (s_iLevel == 2)
	{
		vp[j:j+s_iN-2] = 1.0;
		j += s_iN-1;
	} else
	if (s_iLevel == 3)
	{
		vp[j:j+s_iN-2] = -0.5; // extra log standard deviation level terms
		j += s_iN-1;
	} else
		if (s_iLevel == 4)
	{
		vp[j:j+s_iN-2] = -0.5; // extra log standard deviation level terms
		vp[j+s_iN-2:j+s_iN-2+(s_iN*(s_iN-1))/2] = 1.0;
		j += (s_iN*(s_iN-1))/2 + (s_iN-1);
	}
	//
	if (s_iSlope == 2)
	{
		vp[j:j+s_iN-2] = 1.0;
		j += s_iN-1;
	}
	//
	if (s_iSeaso == 2)
	{
		vp[j:j+s_iN-2] = 1.0;
		j += s_iN-1;
	}
	
	s_vPar = vp;
	return vp;
}

Likelihood1(const vPar, const pdLik, const pvSco, const pmHes)
{
	decl mp, mo, ms;
	[mp,mo,ms] = SetStsmModel1(s_iN, vPar, FALSE);
	SsfLik(pdLik, &s_dVar, s_mY, mp, mo, ms);
	pdLik[0] /= s_iT;         
	return 1;
}
Likelihood1Ex(const vPar, const pdLik, const pvSco, const pmHes)
{
	decl mp, mo, ms;
	[mp,mo,ms] = SetStsmModel1(s_iN, vPar, FALSE);
	decl cst = columns(mp);
	decl csy = rows(mp);
	decl cy = csy - cst;

	// fast likelihood, eq by eq
	if (max(fabs(vPar)) > 500)
    {
        println(" Parameters too large:", vPar);   
    	return 0;
    }
    if (max(fabs(vPar)) <= 500)
    {	SsfLikEx(pdLik, &s_dVar, s_mY, mp, mo, ms);
		pdLik[0] /= s_iT;         
		return 1;
	}
}
SsfPrediction1(const vPar)
{
	decl mp, mo, ms;
	[mp,mo,ms] = SetStsmModel1(s_iN, vPar, FALSE);
	decl cst = columns(mp);
	decl csy = rows(mp);
	decl cy = csy - cst;

	// SsfMomentEstEx, eq by eq
	decl mpred;
	SsfMomentEst(ST_PRED, &mpred, s_mY, mp,mo,ms); // get predictions
	return {mpred[cst:csy-1][],mpred[cst+csy:csy+csy-1][], cst}; // return matrix of predictions Y and their variances and the number of states
}
SsfComponents1(const vPar)
{
	decl mp, mo, ms;
	[mp,mo,ms] = SetStsmModel1(s_iN, vPar, FALSE);
	decl cst = columns(mp);
	decl csy = rows(mp);
	decl cy = csy - cst;
	//println("cy", cy); println("cst", cst); println("csy", csy);

	// SsfMomentEstEx, eq by eq
	decl msmo, mcmp=<>, mcmpvar=<>;
	SsfMomentEstEx(ST_SMO, &msmo, s_mY, mp,mo,ms); // get predictions
	//println("msmo", msmo);
	mcmp = msmo[:cy-1][]; mcmpvar = msmo[csy:csy+cy-1][]; // level estimate and var
	//println("mcmpvar", mcmpvar);
	mcmp |= msmo[cy:cy+cy-1][]; mcmpvar |= msmo[csy+cy:csy+cy+cy-1][]; // slope estimate and var
	mcmp |= mp[cst:][cy+cy:cy+cy+(11*cy)-1] * msmo[cy+cy:cy+cy+(11*cy)-1][]; // seasonal estimate
	mcmp |= msmo[cst-cy:cst-1][]; // irregular estimate
	return {mcmp,mcmpvar}; // return matrix with four cmps (lvl,slp,seas,irr) and their variances
}
SsfForecast1(const vPar)
{
	decl mp, mo, ms;
	[mp,mo,ms] = SetStsmModel1(s_iN, vPar, FALSE);
	decl cst = columns(mp);
	decl csy = rows(mp);
	decl cy = csy - cst;

	// SsfMomentEstEx, eq by eq
	decl msmo, mfor=<>, mforvar=<>;
	SsfMomentEstEx(ST_SMO, &msmo, s_mY, mp,mo,ms); // get predictions

	mfor = msmo[cst:csy-1][];
	mforvar = msmo[csy+cst:csy+csy-1][];
	return {mfor,mforvar}; // return matrix with yhat and their variances
}

EstimateModel1(const fprint, const fgraph, const vIniPar)
{
	decl vpar = vIniPar == <> ? InitPar1() : vIniPar;
//	return vpar;
	if (fprint)
		println("cpar = ", sizerc(vpar), ", initial vpar", vpar); 

	MaxControl(-1, 5, 1);
	decl dlik = 0.0;
	decl ir = MaxBFGS(Likelihood1Ex, &vpar, &dlik, 0, TRUE);
	println( sprint("\n", MaxConvergenceMsg(ir),
          " using numerical derivatives",
	      "\nLog-likelihood = ", dlik * s_iT,
		  "; dVar = ", s_dVar));
		  //, "; \nestimated vpar", vpar));

	s_dLik = dlik * s_iT;
	SetStsmModel1(s_iN, vpar, TRUE); // MODEL 1
	s_vPar = vpar;
	return vpar;
}

ResidualsModeli(const fprint, const fgraph, const iMod, const vPar)
{
	decl mpred, mpredvar, mres, mres0, cst;
	[mpred,mpredvar,cst] = SsfPrediction1(vPar);
	mres = (mpredvar .< 100.0) .? (s_mY - mpred) ./ sqrt(mpredvar) .: M_NAN; // standardize prediction errors
	mres0 = (mres .== M_NAN) .? 0.0 .: mres; // replace missings with 0

	/*if(iMod==1 || iMod==4) s_dAIC = double((-2*s_dLik+2*1)/columns(deletec(mres[0][]))); // calculate AIC of model i --> want alleen level UC, dus 1 relevante univariate hyperparameter (sigma^2_level)
	if(iMod==2 || iMod==3 || iMod==5 || iMod==6) s_dAIC = double((-2*s_dLik+2*3)/columns(deletec(mres[0][]))); // calculate AIC of model i --> level, slope en seasonal UC, dus 3 relevante univariate hyperparameters (sigma^2_level, sigma^2_slope, sigma^2_seasonal)
	*/
	s_cst = cst - s_iN;	// subtract s_iN because irregulars were placed in state vector
	println("total diffuse states:", s_cst);
	s_dAIC = double((-2 * s_dLik + 2 * ((s_cst) + sizerc(vPar))) / s_iNTstar);	// normalized with N*T and total number of diffuse states includes deterministic states, so cst-s_iN because irreg was placed
																										// in state vector. Size of vPar will of course be shorter when states are deterministic.
	s_dBIC = double((-2 * s_dLik + log(s_iT) * ((s_cst) + sizerc(vPar))) / s_iNTstar); 
	
	decl i, k, vcusum, bandcorrgram, bandcusum95, bandcusum90, vones = ones(mres0[0][]), vResstat = <>, vLBstat = <>;
	if (fgraph)
	{
		/*for (i=0; i<s_iN; i++)
		{
			DrawTMatrix(s_igr++, s_mY[i][] | ((mpredvar[i][] .< 100.0) .? mpred[i][] .: M_NAN), {s_sNameY[i], "Prediction"}, s_iFirstYear, s_iFirstMonth, 12, 0, 2);
			//DrawZ(sqrt(mpredvar[i][]), "", ZMODE_BAND, 1.96);
		}
		DrawAdjust(ADJ_AREAMATRIX, 2, 2);
		SaveDrawWindow(sprint("Graphs_v2/PredictionsModel", iMod, ".pdf"));
		ShowDrawWindow();
		CloseDrawWindow();*/
		s_igr=0;
		for (i=0; i<s_iN; i++)
		{
			vcusum = cumulate(mres0[i][]');
			decl iTmiss, iTnomiss;
			iTmiss = sumr(mres[i][:columns(mres)*0.75] .== M_NAN);	// only check for NA's at the start of the time series, assuming NA's only at start and end, but not more than 1/4th at the end and not more than 3/4th at start.
			iTnomiss = columns(deletec(mres[i][]));
			bandcorrgram = 1.96/sqrt(iTnomiss);
			bandcusum95 = .NaN*ones(iTmiss, 1) | 0.948*sqrt(iTnomiss) + 2*0.948*cumulate(vones[:iTnomiss-1]')/sqrt(iTnomiss) | .NaN*ones(s_iT-(iTnomiss+iTmiss),1); // critical values obtained from Harvey (1989) p. 257
			bandcusum90 = .NaN*ones(iTmiss, 1) | 0.850*sqrt(iTnomiss) + 2*0.850*cumulate(vones[:iTnomiss-1]')/sqrt(iTnomiss) | .NaN*ones(s_iT-(iTnomiss+iTmiss),1);
			
			DrawTMatrix(s_igr, mres[i][], sprint("Residual ", s_sNameY[i]), s_iFirstYear, s_iFirstMonth, 12, 0, 2);
			DrawTMatrix(s_igr,  1.96*vones, "", s_iFirstYear, s_iFirstMonth, 12, 0, 3);
			DrawTMatrix(s_igr, -1.96*vones, "", s_iFirstYear, s_iFirstMonth, 12, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Time (month)", 0, 0, -1, -1, TEXT_XLABEL);
			DrawAdjust(ADJ_MINMAX, -4, 4);

			s_igr++;

			DrawDensity(s_igr, deletec(mres[i][]), "Residual", TRUE, TRUE, TRUE, FALSE, FALSE, 20, 2);
			DrawAdjust(ADJ_MINMAX, 0, 0.6);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr++, 0, 200, 0, 1);

			DrawCorrelogram(s_igr, deletec(mres[i][]), "Residual", 40);
			DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			DrawMatrix(s_igr, bandcorrgram*vones[:39], "", 1, 1, 0, 3);
			DrawMatrix(s_igr, -bandcorrgram*vones[:39], "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;
			
			DrawTMatrix(s_igr, vcusum',  sprint("Cusum ",  s_sNameY[i]), s_iFirstYear, s_iFirstMonth, 12, 0, 2);
			DrawTMatrix(s_igr, bandcusum95', "", s_iFirstYear, s_iFirstMonth, 12, 0, 3);
			DrawTMatrix(s_igr, -bandcusum95', "", s_iFirstYear, s_iFirstMonth, 12, 0, 3);
			DrawTMatrix(s_igr, bandcusum90', "", s_iFirstYear, s_iFirstMonth, 12, 0);
			DrawAdjust(ADJ_COLOR, 3, 5);
			DrawTMatrix(s_igr, -bandcusum90', "", s_iFirstYear, s_iFirstMonth, 12, 0);
			DrawAdjust(ADJ_COLOR, 3, 5);
			DrawAdjust(ADJ_MINMAX, -65, 65);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Time (month)", 0, 0, -1, -1, TEXT_XLABEL);
			//DrawAdjust(ADJ_MINMAX, -30, 30);
			s_igr++;
		}
		DrawAdjust(ADJ_AREAMATRIX, 4, 4);
		SaveDrawWindow(sprint("Graphs_v2/ResidualsModel", iMod,".pdf"));
		ShowDrawWindow();
		CloseDrawWindow();
		
		s_igr=0;
		for (i=0; i<s_iN; i++)
		{
			vcusum = cumulate(mres0[i][]');
			decl iTmiss, iTnomiss;
			iTmiss = sumr(mres[i][:columns(mres)*0.75] .== M_NAN);	// only check for NA's at the start of the time series, assuming NA's only at start and end, but not more than 1/4th at the end and not more than 3/4th at start.
			iTnomiss = columns(deletec(mres[i][]));
			bandcusum95 = .NaN*ones(iTmiss, 1) | 0.948*sqrt(iTnomiss) + 2*0.948*cumulate(vones[:iTnomiss-1]')/sqrt(iTnomiss) | .NaN*ones(s_iT-(iTnomiss+iTmiss),1); // critical values obtained from Harvey (1989) p. 257
			bandcusum90 = .NaN*ones(iTmiss, 1) | 0.850*sqrt(iTnomiss) + 2*0.850*cumulate(vones[:iTnomiss-1]')/sqrt(iTnomiss) | .NaN*ones(s_iT-(iTnomiss+iTmiss),1);
			
			DrawTMatrix(s_igr, mres[i][], sprint("Residual ", s_sNameY[i]), s_iFirstYear, s_iFirstMonth, 12, 0, 2);
			DrawTMatrix(s_igr,  1.96*vones, "", s_iFirstYear, s_iFirstMonth, 12, 0, 3);
			DrawTMatrix(s_igr, -1.96*vones, "", s_iFirstYear, s_iFirstMonth, 12, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Time (month)", 0, 0, -1, -1, TEXT_XLABEL);
			DrawAdjust(ADJ_MINMAX, -4, 4);

			s_igr++;

			DrawDensity(s_igr, deletec(mres[i][]), "Residual", TRUE, TRUE, TRUE, FALSE, FALSE, 20, 2);
			DrawAdjust(ADJ_MINMAX, 0, 0.6);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr++, 0, 200, 0, 1);

			DrawQQ(s_igr, deletec(mres[i][]), sprint("Residual ", s_sNameY[i]), QQ_N, 0, 0);
			//DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			//if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;
			
			DrawTMatrix(s_igr, vcusum',  sprint("Cusum ",  s_sNameY[i]), s_iFirstYear, s_iFirstMonth, 12, 0, 2);
			DrawTMatrix(s_igr, bandcusum95', "", s_iFirstYear, s_iFirstMonth, 12, 0, 3);
			DrawTMatrix(s_igr, -bandcusum95', "", s_iFirstYear, s_iFirstMonth, 12, 0, 3);
			DrawTMatrix(s_igr, bandcusum90', "", s_iFirstYear, s_iFirstMonth, 12, 0);
			DrawAdjust(ADJ_COLOR, 3, 5);
			DrawTMatrix(s_igr, -bandcusum90', "", s_iFirstYear, s_iFirstMonth, 12, 0);
			DrawAdjust(ADJ_COLOR, 3, 5);
			DrawAdjust(ADJ_MINMAX, -65, 65);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Time (month)", 0, 0, -1, -1, TEXT_XLABEL);
			//DrawAdjust(ADJ_MINMAX, -30, 30);
			s_igr++;
		}
		DrawAdjust(ADJ_AREAMATRIX, 4, 4);
		SaveDrawWindow(sprint("Graphs_v2/ResidualsModel", iMod,"_v2a.pdf"));
		ShowDrawWindow();
		CloseDrawWindow();

		s_igr=0;
		for (i=0; i<s_iN; i++)
		{
			decl iTmiss, iTnomiss;
			iTmiss = sumr(mres[i][:columns(mres)*0.75] .== M_NAN);	// only check for NA's at the start of the time series, assuming NA's only at start and end, but not more than 1/4th at the end and not more than 3/4th at start.
			iTnomiss = columns(deletec(mres[i][]));
			bandcorrgram = 1.96/sqrt(iTnomiss);
			
			DrawTMatrix(s_igr, sqr(mres[i][]), sprint("Squared residual ", s_sNameY[i]), s_iFirstYear, s_iFirstMonth, 12, 0, 2);
			//DrawTMatrix(s_igr,  1.96*vones, "", s_iFirstYear, s_iFirstMonth, 12, 0, 3);
			//DrawTMatrix(s_igr, -1.96*vones, "", s_iFirstYear, s_iFirstMonth, 12, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Time (month)", 0, 0, -1, -1, TEXT_XLABEL);
			DrawAdjust(ADJ_MINMAX, -4, 4);

			s_igr++;
			DrawCorrelogram(s_igr, deletec(mres[i][]), "Residual", 40);
			DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			DrawMatrix(s_igr, bandcorrgram*vones[:39], "", 1, 1, 0, 3);
			DrawMatrix(s_igr, -bandcorrgram*vones[:39], "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;

			DrawCorrelogram(s_igr, fabs(deletec(mres[i][])), "Absolute residual", 40);
			//DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			//DrawMatrix(s_igr, bandcorrgram*vones[:39], "", 1, 1, 0, 3);
			//DrawMatrix(s_igr, -bandcorrgram*vones[:39], "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;
			
			DrawCorrelogram(s_igr, sqr(deletec(mres[i][])), "Squared residual", 40);
			//DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			//DrawMatrix(s_igr, bandcorrgram*vones[:39], "", 1, 1, 0, 3);
			//DrawMatrix(s_igr, -bandcorrgram*vones[:39], "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;
		}
		DrawAdjust(ADJ_AREAMATRIX, 4, 4);
		SaveDrawWindow(sprint("Graphs_v2/ResidualsModel", iMod,"_v2b.pdf"));
		ShowDrawWindow();
		CloseDrawWindow();
	}
	if (fprint)
	{
		println(sprint("MODEL ", iMod, " DIAGNOSTICS"));
		println(" - loglikehood value = ", s_dLik);
		println(" - sigma2 =            ", s_dVar);
		println(" - AIC =               ", s_dAIC);
		// diagnostics of target series one-step ahead forecast errors
		for (i=0; i<s_iN; i++)
		{
			decl  vresnomiss, iTnomiss, vacf, vacfsq = <>, j, k, dLBstat, dHsked, dm1, dm2, dm3, dm4, dSkew, dzSkew, dKurt, dzKurt, dNorm; // hoe zit het met de tijdsdimensie die verschilt tussen de vier reeksen?
			vresnomiss = deletec(mres[i][]);
			iTnomiss = columns(vresnomiss);
			vacf = acf(vresnomiss', 40); // hier deletec?
			//println("vacf", vacf);
			for (j=0; j<40; j++) {vacfsq |= vacf[j]^2/(iTnomiss-j);}
			//println("vacfsq", vacfsq);
			dLBstat = double(iTnomiss*(iTnomiss+2)*sumc(vacfsq[1:]));	// chi-squared distributed with df: number of lags - (sizerc(vpar)/N) - 1, Harvey (1989) p.259
			println("h for heteroscedasticity test df: ", floor(iTnomiss/3));
			println("omvang mres deletec: ", columns(vresnomiss));
			dHsked = double(sumr(vresnomiss[iTnomiss-1-floor(iTnomiss/3):iTnomiss-1].^2)/sumr(vresnomiss[0:floor(iTnomiss/3)-1].^2)); // F_{h,h} distributed and t starts after initialization period for KF, depends on number of diffuse states per N
			dm1 = double(meanr(vresnomiss));
			dm2 = double(varr(vresnomiss));
			dm3 = double(meanr((vresnomiss-dm1).^3));
			dm4 = double(meanr((vresnomiss-dm1).^4));
			dSkew = double(dm3/sqrt(dm2^3));
			dzSkew = double(dSkew/sqrt(6/iTnomiss));
			println("z-score for skewness: ", dzSkew);
			dKurt = double(dm4/dm2^2);
			dzKurt = double((dKurt-3)/sqrt(24/iTnomiss));
			println("z-score for kurtosis: ", dzKurt);
			dNorm = double(iTnomiss*((dSkew^2/6)+(dKurt-3)^2/24));	// chi-squared distributed with df: 2
			vResstat |= dLBstat | dHsked | dSkew | dKurt | dNorm;
			println("vResstat", vResstat);

			vLBstat |= dLBstat;
			decl vi = <60, 80, 100, 120>;
			foreach (k in vi)
			{
				vacf = acf(vresnomiss', k);
				vacfsq = <>;
				for (j=0; j<k; j++) {vacfsq |= vacf[j]^2/(iTnomiss-j);}
				dLBstat = double(iTnomiss*(iTnomiss+2)*sumc(vacfsq[1:]));
				vLBstat |= dLBstat;
			}
		}		
	}

	return {vResstat, vLBstat};
}
ComponentsModel1(const fprint, const fgraph, const iMod, const vPar)
{
	decl mcmp, mcmpvar, mCI_L, mCI_U;
	[mcmp,mcmpvar] = SsfComponents1(vPar);
	mCI_L = mcmp[:s_iN+s_iN-1][] - 1.96*sqrt(mcmpvar);
	mCI_U = mcmp[:s_iN+s_iN-1][] + 1.96*sqrt(mcmpvar);
	mCI_L[<0,3>][:35] = M_NAN; mCI_U[<0,3>][:35] = M_NAN; // put bounds before 1990 to missing
	mCI_L[<2,4,6,7>][:120] = M_NAN; mCI_U[<2,4,6,7>][:120] = M_NAN; // put bounds before 1997 to missing
	
	
	decl i;
	if (fgraph)
	{ 
		for (i=0; i<s_iN; i++)
		{
			DrawTMatrix(s_igr, s_mY[i][] | mcmp[i][],  {"Series", sprint("Trend ",  s_sNameY[i])}, s_iFirstYear, s_iFirstMonth, 12, 0, 2);
			DrawTMatrix(s_igr, mCI_U[i][], "", s_iFirstYear, s_iFirstMonth, 12, 0, 4);
			DrawTMatrix(s_igr, mCI_L[i][], "", s_iFirstYear, s_iFirstMonth, 12, 0, 4);
			//DrawAxis(s_igr, 0, 3, 5.5, 4, 1, 1, 0, TRUE);
			//DrawAxis(s_igr, 1, s_iFirstYear, 2022, 2020, 20, 20, 12, TRUE);
			DrawAdjust(ADJ_MINMAX, 3, 5.6);
			DrawLegend(s_igr, 50, -130, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr++, 0, 220, 0, 1);
			DrawTMatrix(s_igr, mcmp[s_iN+i][],  sprint("Slope ",  s_sNameY[i]), s_iFirstYear, s_iFirstMonth, 12, 0, 3);
			DrawTMatrix(s_igr, mCI_U[s_iN+i][], "", s_iFirstYear, s_iFirstMonth, 12, 0, 4);
			DrawTMatrix(s_igr, mCI_L[s_iN+i][], "", s_iFirstYear, s_iFirstMonth, 12, 0, 4);
			DrawAdjust(ADJ_MINMAX, -0.00625, 0.0025);
			DrawLegend(s_igr, 50, -130, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr++, 0, 220, 0, 1);
			DrawTMatrix(s_igr, mcmp[s_iN+s_iN+i][],  sprint("Seasonal ",  s_sNameY[i]), s_iFirstYear, s_iFirstMonth, 12, 0, 4);
			DrawAdjust(ADJ_MINMAX, -0.3, 0.2);
			DrawLegend(s_igr, 50, -130, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr++, 0, 220, 0, 1);
			DrawTMatrix(s_igr, mcmp[s_iN+s_iN+s_iN+i][],  sprint("Irregular ",  s_sNameY[i]), s_iFirstYear, s_iFirstMonth, 12, 0, 5);
			DrawAdjust(ADJ_MINMAX, -0.52, 0.5);
			DrawLegend(s_igr, 50, -130, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr++, 0, 220, 0, 1);
		}
		DrawAdjust(ADJ_AREAMATRIX, 4, 4);
		SaveDrawWindow(sprint("Graphs_v2/ComponentsModel", iMod, ".pdf"));
		ShowDrawWindow();
		CloseDrawWindow();
	}
	if (fprint)
	{
		for (i=0; i<s_iN; i++)
		{
			println("FOR VARIABLE ", s_sNameY[i]);
			println("Level  Estimate at 1996M1 and 2020M12 = ", mcmp[i][0], "   ", mcmp[i][s_iT-1]);
			println("Signal Estimate at 1996M1 and 2020M12 = ", mcmp[i][0] + mcmp[s_iN+s_iN+i][0], "   ", mcmp[i][0] + mcmp[s_iN+s_iN+i][s_iT-1]);
			println("sqrt(mcmpvar) trend", sqrt(mcmpvar[i][]));
		}
	}
}
printlatextablediagn(const file, const stitle, const asformat, const m)
{
	decl i, j, ncol = columns(m), nrow = rows(m);
	decl s;
	s = sprint("\\","begin{table}[ht]\n");
	s = sprint(s, "\\","centering\n");
	s = sprint(s, "\\","caption{", stitle, "}\n");
	s = sprint(s, "\\","begin{tabular}{l");
	for (i=0; i<ncol; i++) s = sprint(s, "c");
	s = sprint(s, "}\n");
	s = sprint(s, "\\","hline\n");			  
	for (j=0; j<nrow; j++)
	{
		for (i=0; i<ncol; i++)
			s = (sizeof(asformat)>0) ? sprint(s, " & ", asformat[j], m[j][i]) : sprint(s, " & ", m[j][i]);
		s = (j==0) ? sprint(s, " ", "\\", "\\", "\n", "\\", "hline\n") : sprint(s, " ", "\\", "\\", "\n");
	}
	s = sprint(s, "\\","hline\n");
	s = sprint(s, "\\","end{tabular}\n");
	s = sprint(s, "\\","end{table}\n");

	println(s);
	//if (isfile(file))
		//fprintln(file, s);
}

FirstStudy()
{	
	decl i, iModel=0, mDiagn = <>, vRstats, vLBstats, mLBstats = <>, mPar = zeros(22, 6), mParalt = zeros(22,6);
	println("Variance matrix options for Level/Irregular:  0 == zeros, 1 == ones, 2 == one factor, 3 == diagonal, 4 == full rank");
	println("Variance matrix options for Slope/Seasonal:  0 == zeros, 1 == ones, 2 == one factor");

	for (s_iLevel=1; s_iLevel<3; s_iLevel++)
	{
		for (s_iSlope=0; s_iSlope<3; s_iSlope++)
		{
			s_iSeaso = s_iSlope;
			s_iIrreg = 4;
			//
			iModel++;
			println(sprint("MODEL ", iModel));
			println("Variance matrix settings in this run:");
			println("s_iLevel: ", s_iLevel); println("s_iSlope: ", s_iSlope); println("s_iSeaso: ", s_iSeaso); println("s_iIrreg: ", s_iIrreg);
			// fit proposed model to data
			EstimateModel1(TRUE, FALSE, <>);
			// save vPar of each model setting
			mPar[][iModel-1] = s_vPar;
			if (iModel==1) mParalt[][0] = exp(2.0*s_vPar[0]) | ones(3,1) | zeros(4,1) | zeros(4,1) | s_vVarCorrir;
			if (iModel==2) mParalt[][1] = exp(2.0*s_vPar[0]) | ones(3,1) | exp(2.0*s_vPar[1]) | ones(3,1) | exp(2.0*s_vPar[2]) | ones(3,1) | s_vVarCorrir;
			if (iModel==3) mParalt[][2] = exp(2.0*s_vPar[0]) | ones(3,1) | exp(2.0*s_vPar[1]) | s_vPar[13:15] | exp(2.0*s_vPar[2]) | s_vPar[16:18] | s_vVarCorrir;
			if (iModel==4) mParalt[][3] = exp(2.0*s_vPar[0]) | s_vPar[11:13] | zeros(4,1) | zeros(4,1) | s_vVarCorrir;
			if (iModel==5) mParalt[][4] = exp(2.0*s_vPar[0]) | s_vPar[13:15] | exp(2.0*s_vPar[1]) | ones(3,1) | exp(2.0*s_vPar[2]) | ones(3,1) | s_vVarCorrir;
			if (iModel==6) mParalt[][5] = exp(2.0*s_vPar[0]) | s_vPar[13:15] | exp(2.0*s_vPar[1]) | s_vPar[16:18] | exp(2.0*s_vPar[2]) | s_vPar[19:21] | s_vVarCorrir;
			// residual diagnostic checks
			s_igr = 0;
			[vRstats, vLBstats] = ResidualsModeli(TRUE, TRUE, iModel, s_vPar);
			mDiagn ~= s_dLik |  s_dAIC | s_dBIC | s_iNTstar | sizerc(s_vPar) | s_cst | vRstats;
			mLBstats ~= vLBstats;
			// check smoothed states
			s_igr = 0;
			ComponentsModel1(TRUE, TRUE, iModel, s_vPar);
		}
	}
	s_cMod = iModel;
	decl asformat = new array[rows(mDiagn)+1];
	asformat[0] = sprint("%5d");
	for (i=1; i<=rows(mDiagn); i++) asformat[i] = sprint("%7.3f");
	println("Model diagnostics", mDiagn);
	printlatextablediagn(<>, "Model diagnostics", asformat, range(1,s_cMod) | mDiagn);

	println("mPar", mPar);
	println("mPar alternative", mParalt);
	decl asformatPar = new array[rows(mParalt)+1];
	asformatPar[0] = sprint("%5d");
	for (i=1; i<=rows(mParalt); i++) asformatPar[i] = sprint("%8.4f");
	printlatextablediagn(<>, "Parameter estimates", asformatPar, range(1,s_cMod) | mPar);
	printlatextablediagn(<>, "Parameter estimates", asformatPar, range(1,s_cMod) | mParalt);

	decl asformatLBstats = new array[rows(mLBstats)+1];
	asformatLBstats[0] = sprint("%5d");
	for (i=1; i<=rows(mLBstats); i++) asformatLBstats[i] = sprint("%7.3f");
	printlatextablediagn(<>, "Ljung-Box test statistics for different lag sizes", asformatLBstats, range(1,s_cMod) | mLBstats);

}


ForecastModel1(const fprint, const fgraph, const fparest, const vParEst, const iscenario, const iyear)
{
	decl mybase = s_mY, iTbase = s_iT;
	decl vyobs, vyfor, verr, vyforse;
	decl mfor, mforvar;
	
	decl ct = ((iyear+1) - s_iFirstYear) * 12; // assuming full year, 12 months, year+1 so you go to the first t
											   // of the next year and you deduct 12 t's from there.
	println("ct", ct);
	s_iT = ct;

	// vyobs are the observations that we want to forecast
	vyobs = s_mY[0][ct-12:ct-1];
	
	s_mY = mybase[][:ct-1];
	if (iscenario == 1)
	{
		s_mY[0][ct-12:ct-1] = M_NAN;
		s_mY[1:][] = M_NAN;
	}
	if (iscenario == 2)
	{
		s_mY[:2][ct-12:ct-1] = M_NAN;
		s_mY[3][ct-24:ct-1] = M_NAN;
	}
	if (iscenario == 3)
	{
		s_mY[:1][ct-12:ct-1] = M_NAN;
		s_mY[2][ct-6:ct-1] = M_NAN;
		s_mY[3][ct-24:ct-1] = M_NAN;
	}
	if (iscenario == 4)
	{
		s_mY[:1][ct-12:ct-1] = M_NAN;
		s_mY[2][ct-3:ct-1] = M_NAN;
		s_mY[3][ct-12:ct-1] = M_NAN;
	}
	if (iscenario == 5)
	{
		s_mY[0][ct-12:ct-1] = M_NAN;
		s_mY[3][ct-12:ct-1] = M_NAN;
	}
	
	// estimate parameter if fparest is TRUE OR if vParEst == <>
	decl vpar = fparest ? EstimateModel1(TRUE, FALSE, vParEst) : vParEst;
	if (vpar == <>) vpar = EstimateModel1(TRUE, FALSE, <>);

	[mfor,mforvar] = SsfForecast1(vpar);
	vyfor = mfor[0][ct-12:ct-1];
	verr = vyobs - vyfor;
	vyforse = sqrt(mforvar[0][ct-12:ct-1]);

	decl i, dmse, dmae, dmape, vvec = zeros(3,1);
	dmse  = vvec[0] = double(meanr(sqr(verr)));
	dmae  = vvec[1] = double(meanr(fabs(verr)));
	dmape = vvec[2] = double(100*meanr(fabs(verr./vyobs)));
	
	if (fgraph)
	{
		DrawTMatrix(s_igr, s_mY[0][],  s_sNameY[0], s_iFirstYear, s_iFirstMonth, 12, 0, 1);
		DrawTMatrix(s_igr, .NaN*ones(1,ct-13)~s_mY[0][ct-13]~mfor[0][ct-12:ct-1],  sprint("Forecast ",  s_sNameY[0]), s_iFirstYear, s_iFirstMonth, 12, 0, 2);
		DrawAdjust(ADJ_COLOR, 2, 5);
		DrawZ(.NaN*ones(1,ct-13)~sqrt(mforvar[0][ct-13:ct-1]), "", ZMODE_BAND, 1.96);
		SaveDrawWindow("Graphs_v2/Outofsample.pdf");
		ShowDrawWindow();
	}
	if (fprint)
	{
		println("FORECAST FOR VARIABLE ", s_sNameY[0]);
		println((vyobs|vyfor|vyobs-vyfor|vyforse|(vyobs-vyfor)./vyforse)');
		println("MSE = ", dmse);
		println("MAE = ", dmae);
		println("MAPE = ", dmape);
	}
	s_mY = mybase; s_iT = iTbase;
	return {vvec, verr, vyobs, vyfor, vyforse};
}

printlatextable(const file, const stitle, const iFirstYear, const iLastYear, const asformat, const m)
{
	decl i, j, years, ncol = columns(m), nrow = rows(m);
	decl s;
	s = sprint("\\","begin{table}[ht]\n");
	s = sprint(s, "\\","centering\n");
	s = sprint(s, "\\","caption{", stitle, "}\n");
	s = sprint(s, "\\","begin{tabular}{l");
	for (i=0; i<ncol; i++) s = sprint(s, "c");
	s = sprint(s, "}\n");
	s = sprint(s, "\\","hline\n");
	years = range(iFirstYear, iLastYear);
	s = sprint(s, " & ", "Model");
	for (i=0; i<columns(years); i++) s=sprint(s, " &", "%5d", years[i]);
	s = sprint(s, " & ", "Total");    
	s = sprint(s, " ", "\\", "\\", "\n");
	s = sprint(s, "\\","hline\n");			  
	for (j=0; j<nrow; j++)
	{
		for (i=0; i<ncol; i++)
			s = (sizeof(asformat)>0) ? sprint(s, " & ", asformat[i], m[j][i]) : sprint(s, " & ", m[j][i]);
		s = sprint(s, " ", "\\", "\\", "\n");
	}
	s = sprint(s, "\\","hline\n");
	s = sprint(s, "\\","end{tabular}\n");
	s = sprint(s, "\\","end{table}\n");

	println(s);
	if (isfile(file))
		fprintln(file, s);
}

FullStudy(const iFirstYearFor, const iLastYearFor, const fprint, const fgraph, const fRealTimeParEst)
{
	decl i, k, tyear, iScenario, asformat = new array[iLastYearFor-iFirstYearFor+3], iMod;
	
	decl il = s_iLevel, // 0 == zeros, 1 == ones, 2 == 1 factor, 3 == diagonal 4 == full rank
		 isl = s_iSlope, // 0,1,2
		 ise = s_iSeaso, // 0,1,2
		 ii = s_iIrreg; // 1,2,3,4
		 
	decl file_mse = fopen("table_mse.tex", "w"), file_mae = fopen("table_mae.tex", "w"), file_mape = fopen("table_mape.tex", "w"), file_pars = fopen("table_pars.tex", "w");

	decl vforstat, verr, vyobs, vyfor, vyforse, mForStat = <>, verrtot = <>, vyobstot = <>,	vyfortot = <>, vyforsetot = <>;
	decl myfor = <>, myforse = <>, mMSE = <>, mMAE = <>, mMAPE = <>, mPar = <>;
	
	for (iScenario=1; iScenario<=5; iScenario++)
	{
		for (s_iLevel=1, iMod=0; s_iLevel<3; s_iLevel++)
		{
			for (s_iSlope=0; s_iSlope<3; s_iSlope++)
			{
				s_iSeaso = s_iSlope;
				s_iIrreg = 4;
				iMod++;
				for (tyear=iFirstYearFor, mForStat = <>, verrtot = <>, vyobstot = <>, vyfortot = <>, vyforsetot = <>, s_vPar = <>; tyear<=iLastYearFor; tyear++)
				{
					println("Year in this run: ", tyear);
					println("Forecast scenario in this run: ", iScenario);
					println("MODEL ", iMod);
					println("Variance matrix settings in this run:");
					println("s_iLevel: ", s_iLevel); println("s_iSlope: ", s_iSlope); println("s_iSeaso: ", s_iSeaso); println("s_iIrreg: ", s_iIrreg);
					[vforstat, verr, vyobs, vyfor, vyforse] = ForecastModel1(FALSE, FALSE, fRealTimeParEst, s_vPar, iScenario, tyear);
					mForStat ~= vforstat; verrtot ~= verr; vyobstot ~= vyobs, vyfortot ~= vyfor; vyforsetot ~= vyforse; 
					// save vPar of best model setting
					if (fRealTimeParEst && iScenario==3 && iMod==4)
					{
						decl vPar = <>;
						vPar = exp(2.0*s_vPar[0]) | s_vPar[11:13] | s_vVarCorrir;
						mPar ~= vPar;
						println("mPar", mPar);						
					}					
				}
				mMSE |= mForStat[0][] ~ double(meanr(sqr(verrtot)));
				mMAE |= mForStat[1][] ~ double(meanr(fabs(verrtot)));
				mMAPE |= mForStat[2][] ~ double(100*meanr(fabs(verrtot./vyobstot)));
				myfor |= vyfortot;
				myforse |= vyforsetot;		
			}
		}
	}
	if (fprint)
	{
		asformat[0] = sprint("%5d");
		for (i=1; i<=(iLastYearFor-iFirstYearFor+2); i++) asformat[i] = sprint("%7.3f");
		printlatextable(file_mse, sprint("Forecasting performance"), iFirstYearFor, iLastYearFor, asformat, ones(5,1) ** range(1, s_cMod)' ~ mMSE);
		printlatextable(file_mae, sprint("Forecasting performance"), iFirstYearFor, iLastYearFor, asformat, ones(5,1) ** range(1, s_cMod)' ~ mMAE);
		printlatextable(file_mape, sprint("Forecasting performance"), iFirstYearFor, iLastYearFor, asformat, ones(5,1) ** range(1, s_cMod)' ~ mMAPE);
		for (i=1; i<=(iLastYearFor-iFirstYearFor+1); i++) asformat[i] = sprint("%8.4f");
		printlatextable(file_pars, sprint("Real-time hyperparameter estimates"), iFirstYearFor, iLastYearFor, asformat[1:], mPar);
	}

	fclose(file_mse); fclose(file_mae); fclose(file_mape); fclose(file_pars);

	if (fgraph)
	{
		decl j;
		// scenario comparison plots
		for (j=0; j<=s_cMod-1; j++)
		{
			DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], iFirstYearFor, s_iFirstMonth, 12, 0, 1);
			for (i=0; i<=4; i++)
			{
				DrawTMatrix(s_igr, myfor[i*s_cMod+j][],  sprint("Forecast in scenario ",  i+1), iFirstYearFor, s_iFirstMonth, 12, 0);
				DrawAdjust(ADJ_COLOR, i+2, 5);
				DrawAdjust(ADJ_MINMAX, 3.2, 4.7);
				DrawText(s_igr, "log(RF)", 0, 0, -1, 300 , TEXT_YLABEL);
				DrawText(s_igr, "Time (month)", 0, 0, -1, 300, TEXT_XLABEL);  				
			}
			SaveDrawWindow(sprint("Graphs_v2/ForecastsSWOV_model", j+1, ".pdf"));
			CloseDrawWindow();
		}
		
		for (j=0; j<=s_cMod-1; j++)
		{
			for (i=0; i<=4; i++)
			{
				DrawTMatrix(s_igr, myforse[i*s_cMod+j][],  sprint("Forecast s.e. in scenario ",  i+1), iFirstYearFor, s_iFirstMonth, 12, 0);
				DrawAdjust(ADJ_COLOR, i+2, 1);
				DrawAdjust(ADJ_MINMAX, 0, 0.25);
			}
			SaveDrawWindow(sprint("Graphs_v2/Forecasts_uncertaintySWOV_model", j+1, ".pdf"));
			CloseDrawWindow();
		}
		// model comparison plots
		for (j=0; j<=4; j++)
		{
			DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], iFirstYearFor, s_iFirstMonth, 12, 0, 1);
			for (i=0; i<=s_cMod-1; i++)
			{
				DrawTMatrix(s_igr, myfor[j*s_cMod+i][],  sprint("Forecast model ",  i+1), iFirstYearFor, s_iFirstMonth, 12, 0);
				DrawAdjust(ADJ_COLOR, i+2, 5);
				DrawAdjust(ADJ_MINMAX, 3.2, 4.7);
				DrawText(s_igr, "log(RF)", 0, 0, -1, 300, TEXT_YLABEL);
				DrawText(s_igr, "Time (month)", 0, 0, -1, 300, TEXT_XLABEL);
			}
			SaveDrawWindow(sprint("Graphs_v2/ForecastsSWOV_scenario", j+1, ".pdf"));
			CloseDrawWindow();
		}
		
		for (j=0; j<=4; j++)
		{
			for (i=0; i<=s_cMod-1; i++)
			{
				DrawTMatrix(s_igr, myforse[j*s_cMod+i][],  sprint("Forecast s.e. model ",  i+1), iFirstYearFor, s_iFirstMonth, 12, 0);
				DrawAdjust(ADJ_COLOR, i+2, 1);
				DrawAdjust(ADJ_MINMAX, 0, 0.25);
			}
			SaveDrawWindow(sprint("Graphs_v2/Forecasts_uncertaintySWOV_scenario", j+1, ".pdf"));
			CloseDrawWindow();
		}

		// combined forecast and uncertainty graphs for winning models in scenario 3
		DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], iFirstYearFor, s_iFirstMonth, 12, 0);
		DrawTMatrix(s_igr, myfor[2*s_cMod+4][], "Forecast model 5", iFirstYearFor, s_iFirstMonth, 12, 0, 2);
		DrawZ(myforse[2*s_cMod+4][], "CI", ZMODE_BAND, 1.96);
		DrawAdjust(ADJ_MINMAX, 3.1, 4.5);
		SaveDrawWindow("Graphs_v2/Forecast_model4_scen3_log.pdf");
		CloseDrawWindow();
		DrawTMatrix(s_igr, exp(vyobstot),  s_sNameY[0], iFirstYearFor, s_iFirstMonth, 12, 0, 1);
		DrawTMatrix(s_igr, exp(myfor[2*s_cMod+4][]), "Forecast model 5", iFirstYearFor, s_iFirstMonth, 12, 0, 2);
		SaveDrawWindow("Graphs_v2/Forecast_model5_scen3.pdf");
		CloseDrawWindow();


		
			
	}

	savemat("Output_v2/myfor.mat", myfor);
	savemat("Output_v2/myforse.mat", myforse);
	s_iLevel = il; s_iSlope = isl; s_iSeaso = ise; s_iIrreg = ii;
}

GetGraphs(const fgraph)
{
	if (fgraph)
	{
	 	DrawTMatrix(0, exp(s_mY[0][]), s_sNameY[0], s_iFirstYear, s_iFirstMonth, 12, 0, 1);
		DrawAdjust(ADJ_MINMAX, 0, 200);
	 	DrawTMatrix(1, exp(s_mY[1][]), s_sNameY[1], s_iFirstYear, s_iFirstMonth, 12, 0, 1);
		DrawAdjust(ADJ_MINMAX, 0, 200);
	 	DrawTMatrix(2, exp(s_mY[2][]), s_sNameY[2], s_iFirstYear, s_iFirstMonth, 12, 0, 1);
		DrawAdjust(ADJ_MINMAX, 0, 200);
	 	DrawTMatrix(3, exp(s_mY[3][]), s_sNameY[3], s_iFirstYear, s_iFirstMonth, 12, 0, 1);
		DrawAdjust(ADJ_MINMAX, 0, 200);
		DrawAdjust(ADJ_AREAMATRIX, 2, 2);
		DrawText(4, "RF", 0, 0, -1, 300, TEXT_YLABEL);
		DrawText(4, "Time (month)", 0, 0, -1, 300, TEXT_XLABEL);
		DrawAxisAuto(4,0,FALSE);
		DrawAxisAuto(4,1,FALSE);
		SaveDrawWindow("Graphs_v2/Series.pdf");
		CloseDrawWindow();
		decl i;
		/*for (i=0; i<=12; i++)
		{
			DrawTMatrix(0, exp(s_mY[0][396:419]), "Final (SN)", 2020, s_iFirstMonth, 12, 0, 2);
			DrawAdjust(ADJ_SYMBOLUSE, ST_LINESYMBOLS);
			DrawAdjust(ADJ_SYMBOL, PL_CIRCLE, 75);
		 	DrawTMatrix(0, exp(s_mY[1][396:431]), "Registered (RWS)", 2020, s_iFirstMonth, 12, 0, 4);
			DrawAdjust(ADJ_SYMBOLUSE, ST_LINESYMBOLS);
			DrawAdjust(ADJ_SYMBOL, PL_CIRCLE, 75);
		 	DrawTMatrix(0, exp(s_mY[2][396:(419+i)]), "Dutch citizens (preliminary)", 2020, s_iFirstMonth, 12, 0, 3);
			DrawAdjust(ADJ_SYMBOLUSE, ST_LINESYMBOLS);
			DrawAdjust(ADJ_SYMBOL, PL_CIRCLE, 75);
		 	DrawTMatrix(0, exp(s_mY[3][396:407]), "Dutch citizens (final)", 2020, s_iFirstMonth, 12, 0, 5);
			DrawAdjust(ADJ_COLOR, 5, 2);
			DrawAdjust(ADJ_SYMBOLUSE, ST_LINESYMBOLS);
			DrawAdjust(ADJ_SYMBOL, PL_CIRCLE, 75);
			DrawAdjust(ADJ_MINMAX, 15, 90);
			DrawAdjust(ADJ_LEGEND, 0, 1, 400);
			DrawText(0, "RF", 0, 0, -1, 400, TEXT_YLABEL);
			DrawText(0, "Time (month)", 0, 0, -1, 400, TEXT_XLABEL);
			SaveDrawWindow(sprint("Graphs_v2/Series_zoomed_", i,".pdf"));
			CloseDrawWindow();
		}*/
		decl bandcorrgram, bandcorrgraml1, bandcorrgraml12, bandcorrgraml1l12, vones;
		s_igr=0;
		for (i=0; i<s_iN; i++)
		{
			decl iTmiss, iTnomiss;
			iTnomiss = columns(deletec(s_mY[i][]));
			iTmiss = s_iT-iTnomiss;
			bandcorrgram = 1.96/sqrt(iTnomiss);
			bandcorrgraml1 = 1.96/sqrt(iTnomiss-1);
			bandcorrgraml12 = 1.96/sqrt(iTnomiss-12);
			bandcorrgraml1l12 = 1.96/sqrt(iTnomiss-13);

			DrawCorrelogram(s_igr, deletec(s_mY[i][]), sprint("log(",s_sNameY[i],")"), 48);
			//DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			DrawMatrix(s_igr, bandcorrgram*ones(1,48), "", 1, 1, 0, 3);
			DrawMatrix(s_igr, -bandcorrgram*ones(1,48), "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			//if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;

			DrawCorrelogram(s_igr, deleter(diff(s_mY[i][]',1))', sprint("(1-L)log(",s_sNameY[i],")"), 48);
			//DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			DrawMatrix(s_igr, bandcorrgraml1*ones(1,48), "", 1, 1, 0, 3);
			DrawMatrix(s_igr, -bandcorrgraml1*ones(1,48), "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			//if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;

			/*DrawCorrelogram(s_igr, deleter(diff(s_mY[i][]',12))', sprint("(1-L^12)log(",s_sNameY[i],")"), 48);
			//DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			DrawMatrix(s_igr, bandcorrgram*ones(1,48), "", 1, 1, 0, 3);
			DrawMatrix(s_igr, -bandcorrgram*ones(1,48), "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			//if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;*/

			DrawCorrelogram(s_igr, deleter(diff(diff(s_mY[i][]',12),1))', sprint("(1-L)(1-L^12)log(",s_sNameY[i],")"), 48);
			//DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			DrawMatrix(s_igr, bandcorrgraml1l12*ones(1,48), "", 1, 1, 0, 3);
			DrawMatrix(s_igr, -bandcorrgraml1l12*ones(1,48), "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			//if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;


			DrawCorrelogram(s_igr, deleter(diff(s_mY[i][]',12))', sprint("(1-L^12)log(",s_sNameY[i],")"), 48);
			//DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			DrawMatrix(s_igr, bandcorrgraml12*ones(1,48), "", 1, 1, 0, 3);
			DrawMatrix(s_igr, -bandcorrgraml12*ones(1,48), "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			//if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;
		
			/*DrawCorrelogram(s_igr, deleter(diff(diff(s_mY[i][]',1),1))', sprint("(1-L)(1-L)log(",s_sNameY[i],")"), 48);
			//DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			DrawMatrix(s_igr, bandcorrgram*ones(1,48), "", 1, 1, 0, 3);
			DrawMatrix(s_igr, -bandcorrgram*ones(1,48), "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			//if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;*/

			
			/*DrawCorrelogram(s_igr, sqr(deletec(mres[i][])), "Squared residual", 40);
			//DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			//DrawMatrix(s_igr, bandcorrgram*vones[:39], "", 1, 1, 0, 3);
			//DrawMatrix(s_igr, -bandcorrgram*vones[:39], "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;*/
		}
		DrawAdjust(ADJ_AREAMATRIX, 4, 4);
		SaveDrawWindow(sprint("Graphs_v2/LogSeriesACFs_v2.pdf"));
		ShowDrawWindow();
		CloseDrawWindow();


	}
}

main()
{
	SsfWarning(FALSE);
	format("%12.6g");

	decl mdata = log(loadmat("Gecorrigeerd_Aangevuld_NewMaandcijfers_swov.csv"))';
	
	s_mY = mdata[<1,0,2,3>][2:435]; 				// select and order variables in y_t
	s_sNameY = {"FINRF", "REGRF", "PRERFDC", "FINRFDC"};

	s_iT = columns(s_mY);   				// time series length
	s_iN = rows(s_mY);						// no of time series variable

	decl iTstar, i;
	for (i=0;i<s_iN;i++)
	{
		iTstar = columns(deletec(s_mY[i][]));
		s_iNTstar += iTstar;						// actual no of observations 
	}

	s_iFirstYear = 1987; s_iFirstMonth = 1; // for newmaandcijfers_swov

	println("Variance matrix options for Level/Irregular:  0 == zeros, 1 == ones, 2 == one factor, 3 == diagonal, 4 == full rank");
	println("Variance matrix options for Slope/Seasonal:  0 == zeros, 1 == ones, 2 == one factor");
	println("Variance matrix settings in this run:");
	println("s_iLevel: ", s_iLevel); println("s_iSlope: ", s_iSlope); println("s_iSeaso: ", s_iSeaso); println("s_iIrreg: ", s_iIrreg);

	println("number of observarions : ", s_iT);
	println("number of time series : ", s_iN);
	//println("s_mY",s_mY);

	GetGraphs(TRUE);

	// in-sample

	// model comparison and residual diagnostic checks based on the full dataset	
	s_igr = 0;
	FirstStudy();
	
	// out-of-sample
	
	// full forecasting study for multiple years and multiple model settings
	s_igr = 0;
	s_cMod=6;
	FullStudy(2011, 2022, TRUE, FALSE, TRUE);
	// last argument: re-estimate vpar for each forecast-year TRUE/FALSE

}
