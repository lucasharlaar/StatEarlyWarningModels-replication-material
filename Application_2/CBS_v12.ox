#include <oxstd.h>
#include <oxdraw.h>
#include <packages/ssfpack/ssfpack_ex.h>
#import <maximize>
#include <oxfloat.h>

static decl s_mStsm,
// variance matrix settings
s_iLevel = 2, // 0 == zeros, 1 == ones, 2 == one factor, 3 == diagonal, 4 == two factor, 5 == full rank
s_iSlope = 2, // 0,1,2,3,4,5
s_iSeaso = 2, // 0,1,2,3,4,5
s_iIrreg = 2; // 0 == zeros, 1 == 1 var for yt, 2 == diagonal 1 var for yt, 1 var for xt's
static decl s_mY, s_iN, s_iT, s_iNTstar, s_vK1, s_vK2, s_vK3, s_vPar, s_dVar, s_dLik, s_cst, s_dAIC, s_dBIC, s_cMod, s_vVarCorrir;
static decl s_igr = 0, s_sNameY, s_iFirstYear, s_iFirstQuarter;


StackMonthObs(const maux)
{
	decl iTm=columns(maux)-1, iTq=iTm/3, i, k;	// remove the first observation of 2023, to keep full years
	s_vK1=s_vK2=s_vK3=zeros(1,iTq);
	s_mY = .NaN*ones(4,iTq);
	for (i=0, k=0; i<iTm-2, k<iTq; i+=3, k++)
	{
		s_mY[1][k] = maux[0][i]; s_mY[2][k] = maux[0][i+1]; s_mY[3][k] = maux[0][i+2];
		s_vK1[k] = maux[1][i]; s_vK2[k] = maux[1][i+1]; s_vK3[k] = maux[1][i+2];
	}
}

GetVarianceMatrixDiagonal(const cy, const vpar, const pipar)
{
	pipar[0] = cy;
	if (vpar == <>)
		return <>;
	else
		return diag(exp(2.0 * vpar[:cy-1]));
}

GetVarianceMatrixDiagonal1par(const cy, const vpar, const pipar)
{
	pipar[0] = 1;
	if (vpar == <>)
		return <>;
	else
		return exp(2.0 * vpar)*unit(cy);
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
	decl i,j,k,l;
	// vpar = variance level, variance slope, variance seasonal, variance matrix irreg
	s_mStsm = <CMP_LEVEL,     0,  0, 0;
		       CMP_SLOPE,     0,  0, 0;
	   		   CMP_SEAS_DUMMY, 0, 4, 0>;	// 4 for quarterly data
	k = 0;
	if (s_iLevel != 0) s_mStsm[0][1] = exp(vpar[k++]); // level
	if (s_iSlope != 0) s_mStsm[1][1] = exp(vpar[k++]); // islope 1 == ones var matrix 
	if (s_iSeaso != 0) s_mStsm[2][1] = exp(vpar[k++]); // iseaso 1 == ones var matrix
	j = k+1; // to keep location of sampling error var parameter at k.
	
    decl mphi, momega, msigma, cst;
	GetSsfStsm(s_mStsm, &mphi, &momega, &msigma); // get state space model
	cst = columns(mphi);
	mphi = mphi ** unit(cy);
	msigma = msigma[:cst-1][:cst-1] ** unit(cy);
	msigma |= 0;
	momega = momega ** ones(cy,cy);
	decl m = zeros(cy,cy);
	i = 0;
	if (s_iIrreg == 1)
	{
		m[0][0] = exp(2.0*vpar[j]);
		i = 1;
	}
	else
	if (s_iIrreg == 2)
	{
		m[0][0] = exp(2.0*vpar[j++]);
		m[1:cy-1][1:cy-1] = GetVarianceMatrixDiagonal1par(cy-1, vpar[j], &i);
	}

	momega[cy*cst:][cy*cst:] = m; // all the way at the end

	m = zeros(cy,cy);
	j += i;
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
			println(" .. loadings level = ", (1 | vpar[j:j+cy-2]));
		}
		m = GetVarianceMatrixChol(cy, 1, vpar[l++] | vpar[j:], &i);
		i--;
	} else
	if (s_iLevel == 3)
	{
		m = GetVarianceMatrixDiagonal(cy, vpar[l++] | vpar[j:], &i);
		i--;
	} else
	if (s_iLevel == 4)
	{
		if (fprint)
		{
			println(" .. variance level 1 = ", exp(2.0* vpar[l]));
			println(" .. loadings level 1 = ", (1 | vpar[j+1:j+cy-1]));
			println(" .. variance level 2 = ", exp(2.0* vpar[j]));
			println(" .. loadings level 2 = ", (1 | vpar[j+cy:j+cy+cy-3]));
		}
		m = GetVarianceMatrixChol(cy, 2, vpar[l++] | vpar[j:], &i);
		i--;
	} else
	if (s_iLevel == 5)
	{
		m = GetVarianceMatrixChol(cy, cy, vpar[l++] | vpar[j:], &i);
		i--;
	}

	decl mxt = zeros(cy*cy+cy*cy,columns(s_mY)), mxtrel = ones(4,columns(s_mY)), s, ct_covid=(2020-s_iFirstYear)*4, exfac=1, vecm;
	if (columns(s_mY)>ct_covid-3) mxtrel[0][ct_covid-3] = 1*exfac;
	if (columns(s_mY)>ct_covid-2) mxtrel[0][ct_covid-2] = 1*exfac;
	if (columns(s_mY)>ct_covid-1) mxtrel[0][ct_covid-1] = 1*exfac;
	if (columns(s_mY)>ct_covid) mxtrel[0][ct_covid] = 1*exfac;
	if (columns(s_mY)>ct_covid+1) mxtrel[0][ct_covid+1] = 1*exfac;
	if (columns(s_mY)>ct_covid+2) mxtrel[0][ct_covid+2] = 1*exfac;
	//
	if (columns(s_mY)>ct_covid-3) mxtrel[1][ct_covid-3] = 1;
	if (columns(s_mY)>ct_covid-2) mxtrel[1][ct_covid-2] = 1;
	if (columns(s_mY)>ct_covid-1) mxtrel[1][ct_covid-1] = 1;
	if (columns(s_mY)>ct_covid) mxtrel[1][ct_covid] = 1;
	if (columns(s_mY)>ct_covid+1) mxtrel[1][ct_covid+1] = 1;
	if (columns(s_mY)>ct_covid+2) mxtrel[1][ct_covid+2] = 1;


	vecm = vec(m);
	for (s=0; s<columns(s_mY); s++) mxt[0][s] = mxtrel[0][s]*vecm[0];
	for (s=0; s<columns(s_mY); s++) mxt[1:cy*cy-1][s] = mxtrel[1][s]*vecm[1:cy*cy-1];

	momega[:cy-1][:cy-1] = m;

	j += i;
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
			println(" .. loadings slope = ", (1 | vpar[j:j+cy-2]));
		}
		m = GetVarianceMatrixChol(cy, 1, vpar[l++] | vpar[j:], &i);
		i--;
	} else
	if (s_iSlope == 3)
	{
		m = GetVarianceMatrixDiagonal(cy, vpar[l++] | vpar[j:], &i);
		i--;
	} else
	if (s_iSlope == 4)
	{
		if (fprint)
		{
			println(" .. variance slope 1 = ", exp(2.0* vpar[l]));
			println(" .. loadings slope 1 = ", (1 | vpar[j+1:j+cy-1]));
			println(" .. variance slope 2 = ", exp(2.0* vpar[j]));
			println(" .. loadings slope 2 = ", (1 | vpar[j+cy:j+cy+cy-3]));
		}
		m = GetVarianceMatrixChol(cy, 2, vpar[l++] | vpar[j:], &i);
		i--;
	} else
	if (s_iSlope == 5)
	{
		m = GetVarianceMatrixChol(cy, cy, vpar[l++] | vpar[j:], &i);
		i--;
	}

	if (columns(s_mY)>ct_covid-3) mxtrel[2][ct_covid-3] = 1*exfac;
	if (columns(s_mY)>ct_covid-2) mxtrel[2][ct_covid-2] = 1*exfac;
	if (columns(s_mY)>ct_covid-1) mxtrel[2][ct_covid-1] = 1*exfac;
	if (columns(s_mY)>ct_covid) mxtrel[2][ct_covid] = 1*exfac;
	if (columns(s_mY)>ct_covid+1) mxtrel[2][ct_covid+1] = 1*exfac;
	if (columns(s_mY)>ct_covid+2) mxtrel[2][ct_covid+2] = 1*exfac;
	//
	if (columns(s_mY)>ct_covid-3) mxtrel[3][ct_covid-3] = 1;
	if (columns(s_mY)>ct_covid-2) mxtrel[3][ct_covid-2] = 1;
	if (columns(s_mY)>ct_covid-1) mxtrel[3][ct_covid-1] = 1;
	if (columns(s_mY)>ct_covid) mxtrel[3][ct_covid] = 1;
	if (columns(s_mY)>ct_covid+1) mxtrel[3][ct_covid+1] = 1;
	if (columns(s_mY)>ct_covid+2) mxtrel[3][ct_covid+2] = 1;


	vecm = vec(m);
	for (s=0; s<columns(s_mY); s++) mxt[cy*cy][s] = mxtrel[2][s]*vecm[0];
	for (s=0; s<columns(s_mY); s++) mxt[cy*cy+1:cy*cy+cy*cy-1][s] = mxtrel[3][s]*vecm[1:cy*cy-1];
	
	momega[cy:cy+cy-1][cy:cy+cy-1] = m;
	
	j += i;
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
			println(" .. loadings seasonal = ", (1 | vpar[j:j+cy-2]));
		}
		m = GetVarianceMatrixChol(cy, 1, vpar[l++] | vpar[j:], &i);
		i--;
	} else
	if (s_iSeaso == 3)
	{
		if (fprint)
		{
			println(" .. variances seasonal = ", exp(2.0* vpar[l]) | exp(2.0*vpar[j:j+cy-2]));
		}
		m = GetVarianceMatrixDiagonal(cy, vpar[l++] | vpar[j:], &i);
		i--;
	} else
	if (s_iSeaso == 4)
	{
		if (fprint)
		{
			println(" .. variance seasonal 1 = ", exp(2.0* vpar[l]));
			println(" .. loadings seasonal 1 = ", (1 | vpar[j+1:j+cy-1]));
			println(" .. variance seasonal 2 = ", exp(2.0* vpar[j]));
			println(" .. loadings seasonal 2 = ", (1 | vpar[j+cy:j+cy+cy-3]));
		}
		m = GetVarianceMatrixChol(cy, 2, vpar[l++] | vpar[j:], &i);
		i--;
	} else
	if (s_iSeaso == 5)
	{
		m = GetVarianceMatrixChol(cy, cy, vpar[l++] | vpar[j:], &i);
		i--;
	}
	momega[cy+cy:cy+cy+cy-1][cy+cy:cy+cy+cy-1] = m;
	

	// placing irregular and sampling error in state vector
	cst = columns(mphi); // bigger now, after kronecker product
	decl csy, msampe, phi=0.59; // phi is the autoregressive parameter of the sampling error, obtained from earlier study.  

	m = mphi[cst:][]; // Z_t
	msampe = zeros(1,cy-1) | unit(cy-1);
	mphi = diagcat(diagcat(mphi[:cst-1][:cst-1],phi*unit(cy-1)),zeros(cy,cy)); // diag(T_t,phi,0)
	mphi |= m ~ msampe ~ unit(cy);
	m = momega[cst:][cst:];	// variance matrix irregular
	momega = diagcat(diagcat(diagcat(momega[:cst-1][:cst-1],GetVarianceMatrixDiagonal1par(cy-1,vpar[k],&i)),m), zeros(cy,cy)); // diag(momega,zero matrix), meas eq has no noise
	msigma = diagcat(diagcat(msigma[:cst-1][:cst-1],unit(cy-1)), m) | 0; // diag(msigma,var_samperror,var_irreg) due to scaling of samp error by var(x_{t,i}), we expect its unconditional variance to be close to 1.
	cst = columns(mphi);
	csy = rows(mphi);
	if (fprint)
	{
		println("s_iLevel: ", s_iLevel); println("s_iSlope: ", s_iSlope); println("s_iSeaso: ", s_iSeaso); println("s_iIrreg: ", s_iIrreg);

		println("vpar", vpar);
		println("mOmega level", momega[:cy-1][:cy-1]);
		println("-- rank mOmega level: ", rank(momega[:cy-1][:cy-1]));
		decl corrmlevel, sdilevel;
		sdilevel = 1 ./ sqrt(diagonal(momega[:cy-1][:cy-1]));
		corrmlevel = diag(sdilevel) * momega[:cy-1][:cy-1] * diag(sdilevel);
		println("mOmega level corr mat", corrmlevel);
		println("mOmega slope", momega[cy:cy+cy-1][cy:cy+cy-1]);
		println("-- rank mOmega slope: ", rank(momega[cy:cy+cy-1][cy:cy+cy-1]));
		decl corrmslope, sdislope;
		sdislope = 1 ./ sqrt(diagonal(momega[cy:cy+cy-1][cy:cy+cy-1]));
		corrmslope = diag(sdislope) * momega[cy:cy+cy-1][cy:cy+cy-1] * diag(sdislope);
		println("mOmega slope corr mat", corrmslope);
		println("mOmega seaso", momega[cy+cy:cy+cy+cy-1][cy+cy:cy+cy+cy-1]);
		println("-- rank mOmega seaso: ", rank(momega[cy+cy:cy+cy+cy-1][cy+cy:cy+cy+cy-1]));
		decl corrmseaso, sdiseaso;
		sdiseaso = 1 ./ sqrt(diagonal(momega[cy+cy:cy+cy+cy-1][cy+cy:cy+cy+cy-1]));
		corrmseaso = diag(sdiseaso) * momega[cy+cy:cy+cy+cy-1][cy+cy:cy+cy+cy-1] * diag(sdiseaso);
		println("mOmega seaso corr mat", corrmseaso);
		println("mOmega samp error", momega[cst-cy-cy+1:cst-cy-1][cst-cy-cy+1:cst-cy-1]);
		println("mOmega irreg", momega[cst-cy:cst-1][cst-cy:cst-1]);
		//println("mPhi", mphi);
		//println("mOmega", momega);
		//println("mSigma", msigma);
	}
	// define more time-varying parameter matrices
	decl mJ_phi, mJ_omega, mdelta, mJ_delta, cJ; 
	mxt |= s_vK1 | s_vK2 | s_vK3;
	cJ = rows(mxt);
	mdelta = zeros(csy,1);
	mJ_delta = -1*ones(csy,1);
	mJ_omega = -1*ones(columns(momega),columns(momega));
	for (i=0; i<cy; i++)
		for (j=0; j<cy; j++)
			mJ_omega[j][i] = i*cy+j;
	for (i=4; i<2*cy; i++)
		for (j=4; j<2*cy; j++)
			mJ_omega[j][i] = i*cy+j-4;			
	mJ_phi = -1*ones(csy,cst);
	mJ_phi[csy-3][cst-cy-3] = cJ-3;
	mJ_phi[csy-2][cst-cy-2] = cJ-2;
	mJ_phi[csy-1][cst-cy-1] = cJ-1;
	//if (fprint) {println("mJPhi", mJ_phi); println("mJOmega", mJ_omega);}
	return {mphi, momega, msigma, mdelta, mJ_phi, mJ_omega, mJ_delta, mxt};
}

InitPar1()
{
	decl i, j, k, l, cp;
	decl m;

	// nrpars for Irregular
	if (s_iIrreg == 1) cp = 1;
	else if (s_iIrreg == 2) cp = 2;
	
	cp += 4; // nrpar for variances of level, slope, seasonal, sampling error
	if (s_iLevel == 0) cp--;
	if (s_iSlope == 0) cp--;
	if (s_iSeaso == 0) cp--;
	j = cp;
	i = 0;
	if (s_iLevel == 2)
	{
		m = GetVarianceMatrixChol(s_iN, 1, <>, &i);
		i--; // because the +=4 for the variances was already added.
	} else
	if (s_iLevel == 3)
	{
		m = GetVarianceMatrixDiagonal(s_iN, <>, &i);
		i--;
	} else
	if (s_iLevel == 4)
	{
		m = GetVarianceMatrixChol(s_iN, 2, <>, &i);
		i--;
	} else
	if (s_iLevel == 5)
	{
		m = GetVarianceMatrixChol(s_iN, 4, <>, &i);
		i--;
	}
	//
	k=0;
	if (s_iSlope == 2)
	{
		m = GetVarianceMatrixChol(s_iN, 1, <>, &k);
		k--; // because the +=4 for the variances was already added.
	} else
	if (s_iSlope == 3)
	{
		m = GetVarianceMatrixDiagonal(s_iN, <>, &k);
		k--;
	} else
	if (s_iSlope == 4)
	{
		m = GetVarianceMatrixChol(s_iN, 2, <>, &k);
		k--;
	} else
	if (s_iSlope == 5)
	{
		m = GetVarianceMatrixChol(s_iN, s_iN, <>, &k);
		k--;
	}
	//
	l=0;
	if (s_iSeaso == 2)
	{
		m = GetVarianceMatrixChol(s_iN, 1, <>, &l);
		l--; // because the +=4 for the variances was already added.
	} else
	if (s_iSeaso == 3)
	{
		m = GetVarianceMatrixDiagonal(s_iN, <>, &l);
		l--;
	} else
	if (s_iSeaso == 4)
	{
		m = GetVarianceMatrixChol(s_iN, 2, <>, &l);
		l--;
	} else
	if (s_iSeaso == 5)
	{
		m = GetVarianceMatrixChol(s_iN, s_iN, <>, &l);
		l--;
	}

	cp += i+k+l;

	decl vp = ones(cp,1);
	i=0;
	if (s_iLevel != 0) vp[i++] = -0.5; // log standard deviation level
	if (s_iSlope != 0) vp[i++] = -1.5; // log standard deviation slope 
	if (s_iSeaso != 0) vp[i++] = -1.0; // log standard deviation seasonal
	vp[i++] = -2; // log standard deviation sampling error
	if (s_iIrreg != 0) vp[i++] = 0;
	if (s_iIrreg == 2) vp[i++] = 0;
	//
	if (s_iLevel == 2)
	{
		j += s_iN-1;
	} else
	if (s_iLevel == 3)
	{
		vp[j:j+s_iN-2] = -0.5; // extra log standard deviation level terms
		j += s_iN-1;
	} else
	if (s_iLevel == 4)
	{
		vp[j] = -0.5; // extra log standard deviation level term
		j += (s_iN*(s_iN-1))/2;
	} else
	if (s_iLevel == 5)
	{
		vp[j:j+s_iN-2] = -0.5; // extra log standard deviation level terms
		j += (s_iN*(s_iN-1))/2 + (s_iN-1);
	}
	//
	if (s_iSlope == 2)
	{
		j += s_iN-1;
	} else
	if (s_iSlope == 3)
	{
		vp[j:j+s_iN-2] = -1.5; // extra log standard deviation slope terms
		j += s_iN-1;
	} else
	if (s_iSlope == 4)
	{
		vp[j] = -1.5; // extra log standard deviation slope term
		j += (s_iN*(s_iN-1))/2;
	} else
	if (s_iSlope == 5)
	{
		vp[j:j+s_iN-2] = -1.5; // extra log standard deviation level terms
		j += (s_iN*(s_iN-1))/2 + (s_iN-1);
	}
	//
	if (s_iSeaso == 2)
	{
		j += s_iN-1;
	} else
	if (s_iSeaso == 3)
	{
		vp[j:j+s_iN-2] = -1.0; // extra log standard deviation seasonal terms
		j += s_iN-1;
	} else
	if (s_iSeaso == 4)
	{
		vp[j] = -1.0; // extra log standard deviation seasonal term
		j += (s_iN*(s_iN-1))/2;
	} else
	if (s_iSeaso == 5)
	{
		vp[j:j+s_iN-2] = -1.0; // extra log standard deviation seasonal terms
		j += (s_iN*(s_iN-1))/2 + (s_iN-1);
	}

	s_vPar = vp;
	return vp;
}

Likelihood1Ex(const vPar, const pdLik, const pvSco, const pmHes)
{
	decl mp, mo, ms, md, mjp, mjo, mjd, mx;
	[mp, mo, ms, md, mjp, mjo, mjd, mx] = SetStsmModel1(s_iN, vPar, FALSE);
	decl cst = columns(mp);
	decl csy = rows(mp);
	decl cy = csy - cst;

	// fast likelihood, eq by eq
	/*if (max(fabs(vPar)) > 3000)
    {
        println(" Parameters too large:", vPar);   
    	return 0;
    }
    if (max(fabs(vPar)) <= 3000)
    {*/
		SsfLikEx(pdLik, &s_dVar, s_mY, mp, mo, ms, md, mjp, mjo, mjd, mx);
		pdLik[0] /= s_iT;         
		return 1;
	//}
}

SsfPrediction1(const vPar)
{
	decl mp, mo, ms, md, mjp, mjo, mjd, mx;
	[mp, mo, ms, md, mjp, mjo, mjd, mx] = SetStsmModel1(s_iN, vPar, FALSE);
	decl cst = columns(mp);
	decl csy = rows(mp);
	decl cy = csy - cst;

	// SsfMomentEstEx, eq by eq
	decl mpred;
	SsfMomentEst(ST_PRED, &mpred, s_mY, mp, mo, ms, md, mjp, mjo, mjd, mx); // get predictions
	return {mpred[cst:csy-1][],mpred[cst+csy:csy+csy-1][], cst}; // return matrix of predictions Y and their variances and the number of states
}

SsfComponents1(const vPar)
{
	decl mp, mo, ms, md, mjp, mjo, mjd, mx;
	[mp, mo, ms, md, mjp, mjo, mjd, mx] = SetStsmModel1(s_iN, vPar, FALSE);
	decl cst = columns(mp);
	decl csy = rows(mp);
	decl cy = csy - cst;

	// SsfMomentEstEx, eq by eq
	decl msmo, mcmp=<>, mcmpvar=<>;
	SsfMomentEstEx(ST_SMO, &msmo, s_mY, mp, mo, ms, md, mjp, mjo, mjd, mx); // get predictions

	mcmp = msmo[:cy-1][]; mcmpvar = msmo[csy:csy+cy-1][]; // level estimate and var
	mcmp |= msmo[cy:cy+cy-1][]; mcmpvar |= msmo[csy+cy:csy+cy+cy-1][]; // slope estimate and var
	mcmp |= mp[cst:][cy+cy:cy+cy+(3*cy)-1] * msmo[cy+cy:cy+cy+(3*cy)-1][]; // seasonal estimate
	decl msey = s_vK1 | s_vK2 | s_vK3;
	mcmp |= msey .* msmo[cy+cy+(3*cy):cy+cy+(3*cy)+cy-2][]; // sampling error estimate
	mcmp |= msmo[cy+cy+(3*cy)+cy-1:cy+cy+(3*cy)+cy+cy-2][]; // irregular estimate

	return {mcmp,mcmpvar}; // return matrix with four cmps (lvl,slp,seas,irr) and their variances
}

SsfForecast1(const vPar)
{
	decl mp, mo, ms, md, mjp, mjo, mjd, mx;
	[mp, mo, ms, md, mjp, mjo, mjd, mx] = SetStsmModel1(s_iN, vPar, FALSE);
	decl cst = columns(mp);
	decl csy = rows(mp);
	decl cy = csy - cst;

	// SsfMomentEstEx, eq by eq
	decl msmo, mfor=<>, mforvar=<>;
	SsfMomentEstEx(ST_SMO, &msmo, s_mY,  mp, mo, ms, md, mjp, mjo, mjd, mx); // get predictions

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

	//s_dAIC = double((-2*s_dLik+2*4)/columns(deletec(mres[0][]))); // calculate AIC --> 4 want UC componenten varianties + sampling error variantie
	s_cst = cst - s_iN - (s_iN-1);	// subtract s_iN because irregulars were placed in state vector and subtract s_iN-1 for the sampling error states who are exactly initialised not diffuse
	println("total diffuse states:", s_cst);
	s_dAIC = double((-2 * s_dLik + 2 * ((s_cst) + sizerc(vPar))) / s_iNTstar);	// normalized with N*T and total number of diffuse states includes deterministic states, so cst-s_iN because irreg was placed
																										// in state vector. Size of vPar will of course be shorter when states are deterministic.
	s_dBIC = double((-2 * s_dLik + log(s_iT) * ((s_cst) + sizerc(vPar))) / s_iNTstar); 

	decl i, k, vcusum, bandcorrgram, bandcusum95, bandcusum90, vones = ones(mres0[0][]), vResstat = <>, vLBstat = <>;
	if (fgraph)
	{
		/*for (i=0; i<s_iN; i++)
		{
			DrawTMatrix(s_igr++, s_mY[i][] | ((mpredvar[i][] .< 100.0) .? mpred[i][] .: M_NAN), {s_sNameY[i], "Prediction"}, s_iFirstYear, s_iFirstQuarter, 4, 0, 2);
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
			//iTnomiss = s_iT - iTmiss;
			iTnomiss = columns(deletec(mres[i][]));
			bandcorrgram = 1.96/sqrt(iTnomiss);
			println("iTnomiss", iTnomiss);
			bandcusum95 = .NaN*ones(iTmiss, 1) | 0.948*sqrt(iTnomiss) + 2*0.948*cumulate(vones[:iTnomiss-1]')/sqrt(iTnomiss) | .NaN*ones(s_iT-(iTnomiss+iTmiss),1); // critical values obtained from Harvey (1989) p. 257
			bandcusum90 = .NaN*ones(iTmiss, 1) | 0.850*sqrt(iTnomiss) + 2*0.850*cumulate(vones[:iTnomiss-1]')/sqrt(iTnomiss) | .NaN*ones(s_iT-(iTnomiss+iTmiss),1);
			
			DrawTMatrix(s_igr, mres[i][], sprint("Residual ", s_sNameY[i]), s_iFirstYear, s_iFirstQuarter, 4, 0, 2);
			DrawTMatrix(s_igr,  1.96*vones, "", s_iFirstYear, s_iFirstQuarter, 4, 0, 3);
			DrawTMatrix(s_igr, -1.96*vones, "", s_iFirstYear, s_iFirstQuarter, 4, 0, 3);
			DrawAdjust(ADJ_MINMAX, -5, 5);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 230, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Time (quarter)", 0, 0, -1, -1, TEXT_XLABEL);
			

			s_igr++;

			DrawDensity(s_igr, deletec(mres[i][]), "Residual", TRUE, TRUE, TRUE, FALSE, FALSE, 20, 2);
			DrawAdjust(ADJ_MINMAX, 0, 0.8);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr++, 0, 230, 0, 1);
			
			DrawCorrelogram(s_igr, deletec(mres[i][]), "Residual", 20);
			DrawMatrix(s_igr, bandcorrgram*vones[:19], "", 1, 1, 0, 3);
			DrawMatrix(s_igr, -bandcorrgram*vones[:19], "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 230, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;
			
			DrawTMatrix(s_igr, vcusum',  sprint("Cusum ",  s_sNameY[i]), s_iFirstYear, s_iFirstQuarter, 4, 0, 2);
			DrawTMatrix(s_igr, bandcusum95', "", s_iFirstYear, s_iFirstQuarter, 4, 0, 3);
			DrawTMatrix(s_igr, -bandcusum95', "", s_iFirstYear, s_iFirstQuarter, 4, 0, 3);
			DrawTMatrix(s_igr, bandcusum90', "", s_iFirstYear, s_iFirstQuarter, 4, 0);
			DrawAdjust(ADJ_COLOR, 3, 5);
			DrawTMatrix(s_igr, -bandcusum90', "", s_iFirstYear, s_iFirstQuarter, 4, 0);
			DrawAdjust(ADJ_COLOR, 3, 5);
			DrawAdjust(ADJ_MINMAX, -30, 30);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 230, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Time (quarter)", 0, 0, -1, -1, TEXT_XLABEL);
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
			//iTnomiss = s_iT - iTmiss;
			iTnomiss = columns(deletec(mres[i][]));
			bandcorrgram = 1.96/sqrt(iTnomiss);
			println("iTnomiss", iTnomiss);
			bandcusum95 = .NaN*ones(iTmiss, 1) | 0.948*sqrt(iTnomiss) + 2*0.948*cumulate(vones[:iTnomiss-1]')/sqrt(iTnomiss) | .NaN*ones(s_iT-(iTnomiss+iTmiss),1); // critical values obtained from Harvey (1989) p. 257
			bandcusum90 = .NaN*ones(iTmiss, 1) | 0.850*sqrt(iTnomiss) + 2*0.850*cumulate(vones[:iTnomiss-1]')/sqrt(iTnomiss) | .NaN*ones(s_iT-(iTnomiss+iTmiss),1);
			
			DrawTMatrix(s_igr, mres[i][], sprint("Residual ", s_sNameY[i]), s_iFirstYear, s_iFirstQuarter, 4, 0, 2);
			DrawTMatrix(s_igr,  1.96*vones, "", s_iFirstYear, s_iFirstQuarter, 4, 0, 3);
			DrawTMatrix(s_igr, -1.96*vones, "", s_iFirstYear, s_iFirstQuarter, 4, 0, 3);
			DrawAdjust(ADJ_MINMAX, -5, 5);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 230, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Time (quarter)", 0, 0, -1, -1, TEXT_XLABEL);
			

			s_igr++;

			DrawDensity(s_igr, deletec(mres[i][]), "Residual", TRUE, TRUE, TRUE, FALSE, FALSE, 20, 2);
			DrawAdjust(ADJ_MINMAX, 0, 0.8);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr++, 0, 230, 0, 1);
			
			DrawQQ(s_igr, deletec(mres[i][]), sprint("Residual ", s_sNameY[i]), QQ_N, 0, 0);
			//DrawAdjust(ADJ_MINMAX, -0.25, 0.25);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 200, 0, 1);
			//if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;
			
			DrawTMatrix(s_igr, vcusum',  sprint("Cusum ",  s_sNameY[i]), s_iFirstYear, s_iFirstQuarter, 4, 0, 2);
			DrawTMatrix(s_igr, bandcusum95', "", s_iFirstYear, s_iFirstQuarter, 4, 0, 3);
			DrawTMatrix(s_igr, -bandcusum95', "", s_iFirstYear, s_iFirstQuarter, 4, 0, 3);
			DrawTMatrix(s_igr, bandcusum90', "", s_iFirstYear, s_iFirstQuarter, 4, 0);
			DrawAdjust(ADJ_COLOR, 3, 5);
			DrawTMatrix(s_igr, -bandcusum90', "", s_iFirstYear, s_iFirstQuarter, 4, 0);
			DrawAdjust(ADJ_COLOR, 3, 5);
			DrawAdjust(ADJ_MINMAX, -30, 30);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 230, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Time (quarter)", 0, 0, -1, -1, TEXT_XLABEL);
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
			//iTnomiss = s_iT - iTmiss;
			iTnomiss = columns(deletec(mres[i][]));
			bandcorrgram = 1.96/sqrt(iTnomiss);
			
			DrawTMatrix(s_igr, sqr(mres[i][]), sprint("Squared residual ", s_sNameY[i]), s_iFirstYear, s_iFirstQuarter, 4, 0, 2);
			//DrawTMatrix(s_igr,  1.96*vones, "", s_iFirstYear, s_iFirstQuarter, 4, 0, 3);
			//DrawTMatrix(s_igr, -1.96*vones, "", s_iFirstYear, s_iFirstQuarter, 4, 0, 3);
			//DrawAdjust(ADJ_MINMAX, -5, 5);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 230, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Time (quarter)", 0, 0, -1, -1, TEXT_XLABEL);
			
			s_igr++;
			
			DrawCorrelogram(s_igr, deletec(mres[i][]), "Residual", 20);
			DrawMatrix(s_igr, bandcorrgram*vones[:19], "", 1, 1, 0, 3);
			DrawMatrix(s_igr, -bandcorrgram*vones[:19], "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 230, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;
			
			DrawCorrelogram(s_igr, fabs(deletec(mres[i][])), "Absolute residual", 20);
			//DrawMatrix(s_igr, bandcorrgram*vones[:19], "", 1, 1, 0, 3);
			//DrawMatrix(s_igr, -bandcorrgram*vones[:19], "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 230, 0, 1);
			if (i==s_iN-1) DrawText(s_igr, "Lag", 0, 0, -1, -1, TEXT_XLABEL);
			s_igr++;

			DrawCorrelogram(s_igr, sqr(deletec(mres[i][])), "Squared residual", 20);
			//DrawMatrix(s_igr, bandcorrgram*vones[:19], "", 1, 1, 0, 3);
			//DrawMatrix(s_igr, -bandcorrgram*vones[:19], "", 1, 1, 0, 3);
			DrawLegend(s_igr, 40, 10, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr, 0, 230, 0, 1);
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
			decl  vresnomiss, iTnomiss, vacf, vacfsq = <>, j, dLBstat, dHsked, dm1, dm2, dm3, dm4, dSkew, dzSkew, dKurt, dzKurt, dNorm; // hoe zit het met de tijdsdimensie die verschilt tussen de vier reeksen?
			vresnomiss = deletec(mres[i][]);
			iTnomiss = columns(vresnomiss);
			vacf = acf(vresnomiss', 20); // hier deletec?
			//println("vacf", vacf);
			for (j=0; j<20; j++) {vacfsq |= vacf[j]^2/(iTnomiss-j);}
			//println("vacfsq", vacfsq);
			dLBstat = double(iTnomiss*(iTnomiss+2)*sumc(vacfsq[1:]));	// chi-squared distributed with df: number of lags - number of fitted hyperparameters, Harvey (1989) p.259
			println("h for heteroscedasticity test df: ", floor(iTnomiss/3));
			println("omvang mres deletec: ", columns(vresnomiss));
			dHsked = double(sumr(vresnomiss[iTnomiss-1-floor(iTnomiss/3):iTnomiss-1].^2)/sumr(vresnomiss[0:floor(iTnomiss/3)-1].^2));
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
			dNorm = double(iTnomiss*((dSkew^2/6)+(dKurt-3)^2/24));
			vResstat |= dLBstat | dHsked | dSkew | dKurt | dNorm;
			println("vResstat", vResstat);

			vLBstat |= dLBstat;
			decl vi = <30, 40, 50, 60>;
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
	decl mcmp, mcmpvar;
	[mcmp,mcmpvar] = SsfComponents1(vPar);
	println("mcmp", mcmp);

	decl i;
	if (fgraph)
	{
		for (i=0; i<s_iN; i++)
		{
			DrawTMatrix(s_igr, s_mY[i][] | mcmp[i][],  {"Series", sprint("Trend ",  s_sNameY[i])}, s_iFirstYear, s_iFirstQuarter, 4, 0, 2);
			DrawZ(sqrt(mcmpvar[i][]), "", ZMODE_BAND, 1.96);
			DrawAdjust(ADJ_MINMAX, 60, 72);
			DrawLegend(s_igr, 50, -150, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr++, 0, 270, 0, 1);
			DrawTMatrix(s_igr, mcmp[s_iN+i][],  sprint("Slope ",  s_sNameY[i]), s_iFirstYear, s_iFirstQuarter, 4, 0, 3);
			DrawZ(sqrt(mcmpvar[s_iN+i][]), "", ZMODE_BAND, 1.96);
			DrawAdjust(ADJ_MINMAX, -0.9, 1.2);
			DrawLegend(s_igr, 50, -150, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr++, 0, 270, 0, 1);
			DrawTMatrix(s_igr, mcmp[s_iN+s_iN+i][],  sprint("Seasonal ",  s_sNameY[i]), s_iFirstYear, s_iFirstQuarter, 4, 0, 4);
			DrawAdjust(ADJ_MINMAX, -0.6, 0.5);
			DrawLegend(s_igr, 50, -150, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr++, 0, 270, 0, 1);
			if (i!=0)
			{
				//DrawTMatrix(s_igr, mcmp[s_iN+s_iN+s_iN+i-1][12:],  sprint("Sampling error ",  s_sNameY[i]), s_iFirstYear+3, s_iFirstQuarter, 4, 0, 5); // start in 2004 om rekening te houden met burnin periode kalman filter
				DrawTMatrix(s_igr, mcmp[s_iN+s_iN+s_iN+i-1][],  sprint("Sampling error ",  s_sNameY[i]), s_iFirstYear, s_iFirstQuarter, 4, 0, 5); // start in 2004 om rekening te houden met burnin periode kalman filter
				DrawAdjust(ADJ_MINMAX, -0.5, 0.5);
				DrawLegend(s_igr, 50, -150, FALSE);
				DrawAdjust(ADJ_LEGEND, s_igr++, 0, 270, 0, 1);
			} else
			{
				s_igr++;
			}
			/*DrawTMatrix(s_igr, mcmp[s_iN+s_iN+s_iN+s_iN-1+i][],  sprint("Irregular ",  s_sNameY[i]), s_iFirstYear, s_iFirstQuarter, 4, 0, 5);
			DrawAdjust(ADJ_MINMAX, -0.5, 0.5);
			DrawLegend(s_igr, 50, -130, FALSE);
			DrawAdjust(ADJ_LEGEND, s_igr++, 0, 200, 0, 1);*/
		}
		DrawAdjust(ADJ_AREAMATRIX, 4, 4);
		SaveDrawWindow(sprint("Graphs_v2/ComponentsModel", iMod, ".pdf"));
		ShowDrawWindow();
		CloseDrawWindow();
		/*s_igr=0;
		for (i=0; i<s_iN; i++)
		{
			DrawTMatrix(s_igr++, sqrt(mcmpvar[i][12:]),  sprint("Smoothed level s.e. ",  s_sNameY[i]), s_iFirstYear+3, s_iFirstQuarter, 4, 0, 2); // start in 2004 om rekening te houden met burnin periode kalman filter
		}
		DrawAdjust(ADJ_AREAMATRIX, 4, 1);
		SaveDrawWindow(sprint("Graphs_v2/TimevarylevelseModel", iMod, ".pdf"));
		ShowDrawWindow();
		CloseDrawWindow();

		s_igr=0;
		for (i=4; i<2*s_iN; i++)
		{
			DrawTMatrix(s_igr++, sqrt(mcmpvar[i][12:]),  sprint("Smoothed slope s.e. ",  s_sNameY[i-4]), s_iFirstYear+3, s_iFirstQuarter, 4, 0, 2); // start in 2004 om rekening te houden met burnin periode kalman filter
		}
		DrawAdjust(ADJ_AREAMATRIX, 4, 1);
		SaveDrawWindow(sprint("Graphs_v2/TimevaryslopeseModel", iMod, ".pdf"));
		ShowDrawWindow();
		CloseDrawWindow();

		
		s_igr=0;
		for (i=0; i<s_iN-1; i++)
		{
			DrawTMatrix(s_igr++, mcmp[s_iN+s_iN+s_iN+i][12:],  sprint("Sampling error ",  s_sNameY[i+1]), s_iFirstYear+3, s_iFirstQuarter, 4, 0, 5); // start in 2004 om rekening te houden met burnin periode kalman filter
		}
		DrawAdjust(ADJ_AREAMATRIX, 3, 1);
		SaveDrawWindow(sprint("Graphs_v2/SamplingErrorModel", iMod, ".pdf"));
		ShowDrawWindow();
		CloseDrawWindow();*/

	}
	if (fprint)
	{
		for (i=0; i<s_iN; i++)
		{
			println("FOR VARIABLE ", s_sNameY[i]);
			println("Level  Estimate at 2001Q1 and 2021Q4 = ", mcmp[i][0], "   ", mcmp[i][s_iT-1]);
			println("Signal Estimate at 2001Q1 and 2021Q4 = ", mcmp[i][0] + mcmp[s_iN+s_iN+i][0], "   ", mcmp[i][0] + mcmp[s_iN+s_iN+i][s_iT-1]);
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

printlatextablevpars(const file, const stitle, const asformat, const m)
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
	decl i, iModel=0, loop = <2,4>, mDiagn = <>, vRstats, vLBstats, mLBstats = <>, mPar = zeros(24, 4), mParalt = zeros(24, 4);
	println("Variance matrix options for Level/Slope/Seasonal:  0 == zeros, 1 == ones, 2 == one factor, 3 == diagonal, 4 == two factor, 5 == full rank");
	println("Variance matrix options for Irregular:  0 == zeros, 1 == 1 var for yt, 2 == diagonal with 1 var for yt and 1 var for all xt's");

	foreach (i in loop)
	{
		s_iLevel = i;
		s_iSlope = i;
		s_iIrreg = 2;
		for (s_iSeaso=2; s_iSeaso<=3; s_iSeaso++)
		{
			iModel++;
			println(sprint("MODEL ", iModel));
			println("Variance matrix settings in this run:");
			println("s_iLevel: ", s_iLevel); println("s_iSlope: ", s_iSlope); println("s_iSeaso: ", s_iSeaso); println("s_iIrreg: ", s_iIrreg);
			// fit proposed model to data
			EstimateModel1(TRUE, TRUE, <>);
			// save vPar of each model setting
			mPar[][iModel-1] = s_vPar;
			if (iModel==1) mParalt[][0] = exp(2.0*s_vPar[0])| 0 | s_vPar[6:8] | zeros(2,1) | exp(2.0*s_vPar[1]) | 0 | s_vPar[9:11] | zeros(2,1) | exp(2.0*s_vPar[2]) | zeros(3,1) |s_vPar[12:14] | exp(2.0*s_vPar[3:5]);
			if (iModel==2) mParalt[][1] = exp(2.0*s_vPar[0])| 0 | s_vPar[6:8] | zeros(2,1) | exp(2.0*s_vPar[1]) | 0 | s_vPar[9:11] | zeros(2,1) | exp(2.0*s_vPar[2]) | exp(2.0*s_vPar[12:14]) | zeros(3,1) | exp(2.0*s_vPar[3:5]);
			if (iModel==3) mParalt[][2] = exp(2.0*s_vPar[0])| exp(2.0*s_vPar[6]) | s_vPar[7:11] | exp(2.0*s_vPar[1]) | exp(2.0*s_vPar[12]) | s_vPar[13:17] | exp(2.0*s_vPar[2]) | zeros(3,1) | s_vPar[18:20] | exp(2.0*s_vPar[3:5]);
			if (iModel==4) mParalt[][3] = exp(2.0*s_vPar[0])| exp(2.0*s_vPar[6]) | s_vPar[7:11] | exp(2.0*s_vPar[1]) | exp(2.0*s_vPar[12]) | s_vPar[13:17] | exp(2.0*s_vPar[2]) | exp(2.0*s_vPar[18:20]) | zeros(3,1) | exp(2.0*s_vPar[3:5]);
			// residual diagnostic checks
			s_igr = 0;
			[vRstats, vLBstats] = ResidualsModeli(TRUE, FALSE, iModel, s_vPar);
			mDiagn ~= s_dLik |  s_dAIC | s_dBIC | s_iNTstar | sizerc(s_vPar) | s_cst | vRstats;
			mLBstats ~= vLBstats;
			// check smoothed states
			s_igr = 0;
			ComponentsModel1(TRUE, FALSE, iModel, s_vPar);
		}
	}
	s_cMod = iModel;
	decl asformat = new array[rows(mDiagn)+1];
	asformat[0] = sprint("%5d");
	for (i=1; i<=rows(mDiagn); i++) asformat[i] = sprint("%7.3f");
	println("Model diagnostics", mDiagn);
	printlatextablediagn(<>, "Model diagnostics", asformat, range(1,iModel) | mDiagn);

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

ForecastModel1(const fprint, const fgraph, const fparest, const vParEst, const iscenario, const iyear, const iquarter)
{
	decl mybase = s_mY, iTbase = s_iT, vx1sebase = s_vK1, vx2sebase = s_vK2, vx3sebase = s_vK3;
	decl yobs, yfor, err, yforse;
	decl mfor, mforvar;
	
	decl ct = (((iyear) - s_iFirstYear) * 4)+iquarter; // assuming full year, 4 quarters
	println("ct", ct);
	s_iT = ct;

	// yobs is the observation that we want to forecast
	yobs = s_mY[0][ct-1];
	
	s_mY = mybase[][:ct-1]; s_vK1 = vx1sebase[:ct-1]; s_vK2 = vx2sebase[:ct-1]; s_vK3 = vx3sebase[:ct-1]; 
	if (iscenario == 1)
	{
		s_mY[0][ct-5:ct-1] = M_NAN;
		s_mY[1:][] = M_NAN;
	}
	if (iscenario == 2)
	{
		s_mY[0][ct-5:ct-1] = M_NAN;
		s_mY[2:][ct-1] = M_NAN;
	}
	if (iscenario == 3)
	{
		s_mY[0][ct-5:ct-1] = M_NAN;
		s_mY[3][ct-1] = M_NAN;
	}
	if (iscenario == 4)
	{
		s_mY[0][ct-5:ct-1] = M_NAN;
	}
	
	// estimate parameter if fparest is TRUE OR if vParEst == <>
	decl vpar = fparest ? EstimateModel1(TRUE, FALSE, vParEst) : vParEst;
	if (fparest == FALSE && iquarter == 1) vpar = EstimateModel1(TRUE, FALSE, vParEst);
	if (vpar == <>) vpar = EstimateModel1(TRUE, FALSE, <>);

	[mfor,mforvar] = SsfForecast1(vpar);
	yfor = mfor[0][ct-1];
	err = yobs - yfor;
	yforse = sqrt(mforvar[0][ct-1]);

	if (fgraph)
	{
		DrawTMatrix(s_igr, s_mY[0][],  s_sNameY[0], s_iFirstYear, s_iFirstQuarter, 4, 0, 1);
		DrawTMatrix(s_igr, .NaN*ones(1,ct-6)~s_mY[0][ct-6]~mfor[0][ct-5:ct-1],  sprint("Forecast ",  s_sNameY[0]), s_iFirstYear, s_iFirstQuarter, 4, 0, 2);
		DrawAdjust(ADJ_COLOR, 2, 5);
		DrawZ(.NaN*ones(1,ct-6)~sqrt(mforvar[0][ct-6:ct-1]), "", ZMODE_BAND, 1.96);
		SaveDrawWindow("Graphs_v2/Outofsample.pdf");
		ShowDrawWindow();
	}
	if (fprint)
	{
		println("FORECAST FOR VARIABLE ", s_sNameY[0]);
		println((yobs|yfor|err|yforse|(yobs-yfor)/yforse)');
	}
	
	s_mY = mybase; s_iT = iTbase; s_vK1 = vx1sebase; s_vK2 = vx2sebase; s_vK3 = vx3sebase;
	return {err, yobs, yfor, yforse};
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
	decl i, tyear, tquarter, iScenario, asformat = new array[iLastYearFor-iFirstYearFor+3], loop = <2,4>, iMod;
	
	decl il = s_iLevel, // 0 == zeros, 1 == ones, 2 == one factor, 3 == diagonal, 4 == two factor, 5 == full rank
		 isl = s_iSlope, // 0,1,2,3,4,5
		 ise = s_iSeaso, // 0,1,2,3,4,5
		 ii = s_iIrreg; // 0 == zeros, 1 == 1 var for yt, 2 == diagonal 1 var for yt, 1 var for xt's
		 
	decl file_mse = fopen("table_mse.tex", "w"), file_mae = fopen("table_mae.tex", "w"), file_mape = fopen("table_mape.tex", "w"), file_pars = fopen("table_pars.tex", "w");

	decl vForStat=zeros(3,1), err, yobs, yfor, yforse, iAIC, verr, vyobs, vyfor, vyforse, v;
	decl mForStat = <>, verrtot = <>, vyobstot = <>, vyfortot = <>, vyforsetot = <>;
	decl myfor = <>, myforse = <>, mMSE = <>, mMAE = <>, mMAPE = <>, mPar = <>;
	
	for (iScenario=1; iScenario<=4; iScenario++)
	{
		iMod=0;
		foreach (i in loop)
		{
			s_iLevel = i;
			s_iSlope = i;
			for (s_iSeaso=2; s_iSeaso<=3; s_iSeaso++)
			{
				//s_iSeaso = i;
				s_iIrreg = 2;
				iMod++;
				for (tyear=iFirstYearFor, mForStat = <>, verrtot = <>, vyobstot = <>, vyfortot = <>, vyforsetot = <>, s_vPar = <>; tyear<=iLastYearFor; tyear++)
				{
					for (tquarter = 1, verr = <>, vyobs = <>, vForStat = <>, vyfor = <>, vyforse = <>; tquarter<=4; tquarter++)	 // wanneer moet s_vPar leeg?
					{
						println("Year and quarter in this run: ", tyear, "Q", tquarter);
						println("Forecast scenario in this run: ", iScenario);
						println("MODEL ", iMod);
						println("Variance matrix settings in this run:");
						println("s_iLevel: ", s_iLevel); println("s_iSlope: ", s_iSlope); println("s_iSeaso: ", s_iSeaso); println("s_iIrreg: ", s_iIrreg);
						[err, yobs, yfor, yforse] = ForecastModel1(FALSE, FALSE, fRealTimeParEst, s_vPar, iScenario, tyear, tquarter);
						verr ~= err; vyobs ~= yobs;	vyfor ~= yfor; vyforse ~= yforse;
					}
					if (tyear != 2022){
					vForStat = double(meanr(sqr(verr))) | double(meanr(fabs(verr))) | double(100*meanr(fabs(verr./vyobs)));
					mForStat ~= vForStat; verrtot ~= verr; vyobstot ~= vyobs;
					}
					vyfortot ~= vyfor; vyforsetot ~= vyforse;
				}
				mMSE |= mForStat[0][] ~ double(meanr(sqr(verrtot)));
				mMAE |= mForStat[1][] ~ double(meanr(fabs(verrtot)));
				mMAPE |= mForStat[2][] ~ double(100*meanr(fabs(verrtot./vyobstot)));
				myfor |= vyfortot;
				myforse |= vyforsetot;
			}
		}
	}

	decl vMSEbase=mMSE[0][], vMAEbase=mMAE[0][], vMAPEbase=mMAPE[0][];
	mMSE ./= vMSEbase;
	mMAE ./= vMAEbase;
	mMAPE ./= vMAPEbase;

	if (fprint)
	{
		asformat[0] = sprint("%5d");
		for (i=1; i<=(iLastYearFor-1-iFirstYearFor+2); i++) asformat[i] = sprint("%7.3f");
		printlatextable(file_mse, sprint("Forecasting performance"), iFirstYearFor, iLastYearFor-1, asformat, ones(4,1) ** range(1, 4)' ~ mMSE);
		printlatextable(file_mae, sprint("Forecasting performance"), iFirstYearFor, iLastYearFor-1, asformat, ones(4,1) ** range(1, 4)' ~ mMAE);
		printlatextable(file_mape, sprint("Forecasting performance"), iFirstYearFor, iLastYearFor-1, asformat, ones(4,1) ** range(1, 4)' ~ mMAPE);
	}

	fclose(file_mse); fclose(file_mae); fclose(file_mape);

	if (fgraph)
	{
		decl j;
		// scenario comparison plots
		for (j=0; j<=s_cMod-1; j++)
		{
			DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], iFirstYearFor, s_iFirstQuarter, 4, 0, 1);
			for (i=0; i<=3; i++)
			{
				DrawTMatrix(s_igr, myfor[i*4+j][],  sprint("Forecast in scenario ",  i+1), iFirstYearFor, s_iFirstQuarter, 4, 0);
				DrawAdjust(ADJ_COLOR, i+2, 5);
			}
			SaveDrawWindow(sprint("Graphs_v2/Forecasts_model", j+1, ".pdf"));
			CloseDrawWindow();
		}
		
		for (j=0; j<=s_cMod-1; j++)
		{
			for (i=0; i<=3; i++)
			{
				DrawTMatrix(s_igr, myforse[i*4+j][],  sprint("Forecast s.e. in scenario ",  i+1), iFirstYearFor, s_iFirstQuarter, 4, 0, i+2);
				DrawAdjust(ADJ_COLOR, i+2, 1);
			}
			SaveDrawWindow(sprint("Graphs_v2/Forecasts_uncertainty_model", j+1, ".pdf"));
			CloseDrawWindow();
		}
		// model comparison plots
		for (j=0; j<=3; j++)
		{
			DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], iFirstYearFor, s_iFirstQuarter, 4, 0, 1);
			for (i=0; i<=s_cMod-1; i++)
			{
				DrawTMatrix(s_igr, myfor[j*s_cMod+i][],  sprint("Forecast model ",  i+1), iFirstYearFor, s_iFirstQuarter, 4, 0);
				DrawAdjust(ADJ_COLOR, i+2, 5);
			}
			SaveDrawWindow(sprint("Graphs_v2/Forecasts_scenario", j+1, ".pdf"));
			CloseDrawWindow();
		}
		
		for (j=0; j<=3; j++)
		{
			for (i=0; i<=s_cMod-1; i++)
			{
				DrawTMatrix(s_igr, myforse[j*s_cMod+i][],  sprint("Forecast s.e. model ",  i+1), iFirstYearFor, s_iFirstQuarter, 4, 0, i+2);
				DrawAdjust(ADJ_COLOR, i+2, 1);
			}
			SaveDrawWindow(sprint("Graphs_v2/Forecasts_uncertainty_scenario", j+1, ".pdf"));
			CloseDrawWindow();
		}

		// combined forecast and uncertainty graph for model 3, scenario 3
		DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], iFirstYearFor, s_iFirstQuarter, 4, 0, 1);
		DrawTMatrix(s_igr, myfor[2*s_cMod+2][], "Forecast in scenario 3", iFirstYearFor, s_iFirstQuarter, 4, 0, 2);
		DrawZ(myforse[2*s_cMod+2][], "CI", ZMODE_BAND, 1.96);
		SaveDrawWindow("Graphs_v2/Forecasts_model3_scen3.pdf");
		CloseDrawWindow();
			
	}

	savemat("Output_v2/myfor.mat", myfor);
	savemat("Output_v2/myforse.mat", myforse);
	s_iLevel = il; s_iSlope = isl; s_iSeaso = ise; s_iIrreg = ii;
}

GetGraphs(const fgraph, const mx)
{
	if (fgraph)
	{
	 	/*DrawTMatrix(0, s_mY[0][], "Register", s_iFirstYear, s_iFirstQuarter, 4, 0, 1);
		DrawAdjust(ADJ_MINMAX, 59.75, 72.5);
	 	DrawTMatrix(1, mx[0][:columns(mx)-2]~M_NAN, "LFS", s_iFirstYear, 1, 12, 0, 1);
		DrawAdjust(ADJ_MINMAX, 59.75, 72.5);
		SaveDrawWindow("Graphs_v2/Series_CBS.pdf");
		ShowDrawWindow();*/
		DrawTMatrix(0, s_mY[0][], "Register Amsterdam (quarterly)", s_iFirstYear, s_iFirstQuarter, 4, 0, 2);
		DrawAdjust(ADJ_SYMBOLUSE, ST_LINESYMBOLS);
		//DrawAdjust(ADJ_SYMBOLUSE, ST_SYMBOLS);
		//DrawAdjust(ADJ_MINMAX, 63, 72);
		DrawAdjust(ADJ_SYMBOL, PL_CIRCLE, 40);
		//DrawAxisAuto(0,1,TRUE,ANCHOR_USER);
		//DrawAxis(0,1);
	 	DrawTMatrix(0, mx[0][:columns(mx)-2]~M_NAN, "LFS Netherlands (monthly)", s_iFirstYear, 1, 12, 0, 3);
		DrawAdjust(ADJ_SYMBOLUSE, ST_LINESYMBOLS);
		DrawAdjust(ADJ_SYMBOL, PL_CIRCLE, 40);
		DrawAdjust(ADJ_MINMAX, 59.75, 72.5);
		//DrawAxisAuto(1,1,TRUE,ANCHOR_MIN);
		SaveDrawWindow("Graphs_v2/Series_CBS.pdf");
		ShowDrawWindow();

		/*DrawTMatrix(0, s_mY[0][56:], "Register Amsterdam (quarterly)", 2015, s_iFirstQuarter, 4, 0, 2);
		DrawAdjust(ADJ_SYMBOLUSE, ST_LINESYMBOLS);
		//DrawAdjust(ADJ_SYMBOLUSE, ST_SYMBOLS);
		//DrawAdjust(ADJ_MINMAX, 63, 72);
		DrawAdjust(ADJ_SYMBOL, PL_CIRCLE, 75);
		//DrawAxisAuto(0,1,TRUE,ANCHOR_USER);
		//DrawAxis(0,1);
	 	DrawTMatrix(0, mx[0][168:columns(mx)-2]~M_NAN, "LFS Netherlands (monthly)", 2015, 1, 12, 0, 3);
		DrawAdjust(ADJ_SYMBOLUSE, ST_LINESYMBOLS);
		DrawAdjust(ADJ_SYMBOL, PL_CIRCLE, 75);
		DrawAdjust(ADJ_MINMAX, 64, 72);
		//DrawAxisAuto(1,1,TRUE,ANCHOR_MIN);
		SaveDrawWindow("Graphs_v2/Series_CBS_zoomed.pdf");
		ShowDrawWindow();*/

		decl iTmiss, iTnomiss, bandcorrgram, bandcorrgraml1, bandcorrgraml4, bandcorrgraml12, bandcorrgraml1l4, bandcorrgraml1l12;
		iTnomiss = columns(deletec(s_mY[0][]));
		iTmiss = s_iT-iTnomiss;
		bandcorrgram = 1.96/sqrt(iTnomiss);
		bandcorrgraml1 = 1.96/sqrt(iTnomiss-1);
		bandcorrgraml4 = 1.96/sqrt(iTnomiss-4);
		bandcorrgraml1l4 = 1.96/sqrt(iTnomiss-5);
		s_igr = 0;
		println("iTnomiss: ", iTnomiss);

		DrawCorrelogram(s_igr, deletec(s_mY[0][]), "Register", 16);
		DrawMatrix(s_igr, bandcorrgram*ones(1,16), "", 1, 1, 0, 3);
		DrawMatrix(s_igr, -bandcorrgram*ones(1,16), "", 1, 1, 0, 3);
		DrawLegend(s_igr, 40, 10, FALSE);
		DrawAdjust(ADJ_LEGEND, s_igr, 200, 0, 1);
		s_igr++;

		DrawCorrelogram(s_igr, deleter(diff(s_mY[0][]',1))', "(1-L)Register", 16);
		DrawMatrix(s_igr, bandcorrgraml1*ones(1,16), "", 1, 1, 0, 3);
		DrawMatrix(s_igr, -bandcorrgraml1*ones(1,16), "", 1, 1, 0, 3);
		DrawLegend(s_igr, 40, 10, FALSE);
		DrawAdjust(ADJ_LEGEND, s_igr, 200, 0, 1);
		s_igr++;

		DrawCorrelogram(s_igr, deleter(diff(diff(s_mY[0][]',4),1))', "(1-L)(1-L^4)Register", 16);
		DrawMatrix(s_igr, bandcorrgraml1l4*ones(1,16), "", 1, 1, 0, 3);
		DrawMatrix(s_igr, -bandcorrgraml1l4*ones(1,16), "", 1, 1, 0, 3);
		DrawLegend(s_igr, 40, 10, FALSE);
		DrawAdjust(ADJ_LEGEND, s_igr, 200, 0, 1);
		s_igr++;

		DrawCorrelogram(s_igr, deleter(diff(s_mY[0][]',4))', "(1-L^4)Register", 16);
		DrawMatrix(s_igr, bandcorrgraml4*ones(1,16), "", 1, 1, 0, 3);
		DrawMatrix(s_igr, -bandcorrgraml4*ones(1,16), "", 1, 1, 0, 3);
		DrawLegend(s_igr, 40, 10, FALSE);
		DrawAdjust(ADJ_LEGEND, s_igr, 200, 0, 1);
		s_igr++;

		iTnomiss = columns(deletec(mx[0][:columns(mx)-2]));
		bandcorrgram = 1.96/sqrt(iTnomiss);
		bandcorrgraml1 = 1.96/sqrt(iTnomiss-1);
		bandcorrgraml12 = 1.96/sqrt(iTnomiss-12);
		bandcorrgraml1l12 = 1.96/sqrt(iTnomiss-13);

		DrawCorrelogram(s_igr, deletec(mx[0][:columns(mx)-2]), "LFS", 48);
		DrawMatrix(s_igr, bandcorrgram*ones(1,48), "", 1, 1, 0, 3);
		DrawMatrix(s_igr, -bandcorrgram*ones(1,48), "", 1, 1, 0, 3);
		DrawLegend(s_igr, 40, 10, FALSE);
		DrawAdjust(ADJ_LEGEND, s_igr, 200, 0, 1);
		s_igr++;

		println("mx: ", deleter(diff(mx[0][:columns(mx)-2]',1))');
		println("iTnomiss: ", iTnomiss);
		DrawCorrelogram(s_igr, deleter(diff(mx[0][:columns(mx)-2]',1))', "(1-L)LFS", 48);
		DrawMatrix(s_igr, bandcorrgraml1*ones(1,48), "", 1, 1, 0, 3);
		DrawMatrix(s_igr, -bandcorrgraml1*ones(1,48), "", 1, 1, 0, 3);
		DrawLegend(s_igr, 40, 10, FALSE);
		DrawAdjust(ADJ_LEGEND, s_igr, 200, 0, 1);
		s_igr++;

		DrawCorrelogram(s_igr, deleter(diff(diff(mx[0][:columns(mx)-2]',12),1))', "(1-L)(1-L^12)LFS", 48);
		DrawMatrix(s_igr, bandcorrgraml1l12*ones(1,48), "", 1, 1, 0, 3);
		DrawMatrix(s_igr, -bandcorrgraml1l12*ones(1,48), "", 1, 1, 0, 3);
		DrawLegend(s_igr, 40, 10, FALSE);
		DrawAdjust(ADJ_LEGEND, s_igr, 200, 0, 1);
		s_igr++;

		DrawCorrelogram(s_igr, deleter(diff(mx[0][:columns(mx)-2]',12))', "(1-L^12)LFS", 48);
		DrawMatrix(s_igr, bandcorrgraml12*ones(1,48), "", 1, 1, 0, 3);
		DrawMatrix(s_igr, -bandcorrgraml12*ones(1,48), "", 1, 1, 0, 3);
		DrawLegend(s_igr, 40, 10, FALSE);
		DrawAdjust(ADJ_LEGEND, s_igr, 200, 0, 1);

		DrawAdjust(ADJ_AREAMATRIX, 2, 4);
		SaveDrawWindow(sprint("Graphs_v2/SeriesACFs_v2.pdf"));
		ShowDrawWindow();
		CloseDrawWindow();


	}
}

main()
{
	SsfWarning(FALSE);
	format("%12.6g");

	decl mx = loadmat("EarlyWarningExampleWerkzameBeroepsBLucas.xlsx")[1:][2:]', vy = loadmat("Amsterdam_Ox_Nap_RegisterKwartaal_Nieuw.xls")';
	StackMonthObs(mx);

	s_mY[0][8:] = vy; 		//register series starts two years later
	s_sNameY = {"Register", "LFS M1", "LFS M2", "LFS M3"};

	s_iT = columns(s_mY);   				// time series length
	s_iN = rows(s_mY);						// no of time series variable

	decl iTstar, i;
	for (i=0;i<s_iN;i++)
	{
		iTstar = columns(deletec(s_mY[i][]));
		s_iNTstar += iTstar;						// actual no of observations 
	}


	println("number of observarions : ", s_iT);
	println("number of time series : ", s_iN);
	
	s_iFirstYear = 2001; s_iFirstQuarter = 1;	// needed for nice graphs

	GetGraphs(TRUE,mx);					// visualize time series

	// in-sample

	// model comparison and residual diagnostic checks based on the full dataset
	s_igr = 0;
	FirstStudy();
	
	// out-of-sample
	
	// full forecasting study for multiple years/quarters
	s_igr = 0;
	s_cMod = 4;
	FullStudy(2011, 2022, TRUE, FALSE, TRUE);
	// last argument: re-estimate vpar for each nowcast-quarter TRUE/FALSE
	
}
