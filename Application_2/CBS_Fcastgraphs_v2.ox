#include <oxstd.oxh>
#include <oxdraw.h>

static decl s_mY, s_iN, s_iT, s_vK1, s_vK2, s_vK3, s_vPar, s_dVar, s_dLik, s_dAIC, s_cMod, s_vVarCorrir;
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

main()
{
	decl mx = loadmat("EarlyWarningExampleWerkzameBeroepsBLucas.xlsx")[1:][2:]', vy = loadmat("Amsterdam_Ox_Nap_RegisterKwartaal_Nieuw.xls")';
	StackMonthObs(mx);

	s_mY[0][8:] = vy; 		//register series starts two years later
	s_sNameY = {"LFP rate register", "LFP rate LFS M1", "LFP rate LFS M2", "LFP rate LFS M3"};

	decl myfor = loadmat("Output_v2/myfor.mat"), myforCovid = loadmat("Output_v2/myforCovid.mat"), myforse = loadmat("Output_v2/myforse.mat"), myforseCovid = loadmat("Output_v2/myforseCovid.mat");
	decl vyobstot = s_mY[0][40:];
	println("vyobstot: ", vyobstot);
	println("myfor: ", myfor);
	println("myforCovid: ", myforCovid);
	println("myforse: ", myforse);
	println("s_mY: ", s_mY);
	/*savemat("vyobstot.zip", vyobstot');
	savemat("myfor.zip", myfor');
	savemat("myforse.zip", myforse');
	savemat("myforCovid.zip", myforCovid');
	savemat("myforCovidse.zip", myforseCovid');*/

	s_iFirstYear = 2001; s_iFirstQuarter = 1;
	s_cMod=4;
	decl i, j;
		// scenario comparison plots
		/*for (j=0; j<=s_cMod-1; j++)
		{
			DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, s_iFirstQuarter, 4, 0, 1);
			for (i=0; i<=3; i++)
			{
				DrawTMatrix(s_igr, myfor[i*4+j][],  sprint("Forecast in scenario ",  i+1), 2011, s_iFirstQuarter, 4, 0);
				DrawAdjust(ADJ_COLOR, i+2, 5);
				DrawAdjust(ADJ_MINMAX, 59.75, 72.25);
			}
			s_igr++;
			DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, s_iFirstQuarter, 4, 0, 1);
			for (i=0; i<=3; i++)
			{
				DrawTMatrix(s_igr, myforCovid[i*4+j][],  sprint("Forecast in scenario ",  i+1), 2011, s_iFirstQuarter, 4, 0);
				DrawAdjust(ADJ_COLOR, i+2, 5);
				DrawAdjust(ADJ_MINMAX, 59.75, 72.25);
			}
			DrawAdjust(ADJ_AREAMATRIX, 2, 1);
			SaveDrawWindow(sprint("Graphs_v2/Forecasts_CBS_model", j+1, ".pdf"));
			CloseDrawWindow();
			s_igr=0;
		}*/
		/*j=2;
		//for (j=0; j<=s_cMod-1; j++)
		//{
			for (i=0; i<=3; i++)
			{
				DrawTMatrix(s_igr, myforse[i*4+j][],  sprint("Forecast s.e. in scenario ",  i+1), 2011, s_iFirstQuarter, 4, 0, i+2);
				DrawAdjust(ADJ_COLOR, i+2, 1);
				DrawAdjust(ADJ_MINMAX, 0.4, 2.8);
			}
			s_igr++;
			for (i=0; i<=3; i++)
			{
				DrawTMatrix(s_igr, myforseCovid[i*4+j][],  sprint("Forecast s.e. in scenario ",  i+1), 2011, s_iFirstQuarter, 4, 0, i+2);
				DrawAdjust(ADJ_COLOR, i+2, 1);
				DrawAdjust(ADJ_MINMAX, 0.4, 2.8);
			}
			s_igr++;
			DrawTMatrix(s_igr, myforse[0][], "Forecast s.e. univariate ", 2011, s_iFirstQuarter, 4, 0, 2);
			for (i=0; i<=s_cMod-1; i++)
			{
				DrawTMatrix(s_igr, myforse[j*s_cMod+i][],  sprint("Forecast s.e. model ",  i+1), 2011, s_iFirstQuarter, 4, 0, i+3);
				DrawAdjust(ADJ_COLOR, i+3, 1);
				DrawAdjust(ADJ_MINMAX, 0.4, 2.8);
			}
			s_igr++;
			DrawTMatrix(s_igr, myforseCovid[0][], "Forecast s.e. univariate ", 2011, s_iFirstQuarter, 4, 0, 2);
			for (i=0; i<=s_cMod-1; i++)
			{
				DrawTMatrix(s_igr, myforseCovid[j*s_cMod+i][],  sprint("Forecast s.e. model ",  i+1), 2011, s_iFirstQuarter, 4, 0, i+3);
				DrawAdjust(ADJ_COLOR, i+3, 1);
				DrawAdjust(ADJ_MINMAX, 0.4, 2.8);
			}
			DrawAdjust(ADJ_AREAMATRIX, 2, 2);
			SaveDrawWindow(sprint("Graphs_v2/Forecasts_uncertainty_CBS.pdf"));
			CloseDrawWindow();*/
			s_igr=0;
		//}
		// model comparison plots
		/*for (j=1; j<=3; j++)
		{
			DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, s_iFirstQuarter, 4, 0, 1);
			DrawTMatrix(s_igr, myfor[0][],  "Forecast univariate", 2011, s_iFirstQuarter, 4, 0, 2);
			for (i=0; i<=s_cMod-1; i++)
			{
				DrawTMatrix(s_igr, myfor[j*s_cMod+i][],  sprint("Forecast model ",  i+1), 2011, s_iFirstQuarter, 4, 0);
				DrawAdjust(ADJ_COLOR, i+3, 5);
				DrawAdjust(ADJ_MINMAX, 62, 71.5);
			}
			s_igr++;
			DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, s_iFirstQuarter, 4, 0, 1);
			DrawTMatrix(s_igr, myforCovid[0][],  "Forecast univariate", 2011, s_iFirstQuarter, 4, 0, 2);
			for (i=0; i<=s_cMod-1; i++)
			{
				DrawTMatrix(s_igr, myforCovid[j*s_cMod+i][],  sprint("Forecast model ",  i+1), 2011, s_iFirstQuarter, 4, 0);
				DrawAdjust(ADJ_COLOR, i+3, 5);
				DrawAdjust(ADJ_MINMAX, 62, 71.5);
			}
			DrawAdjust(ADJ_AREAMATRIX, 2, 1);
			SaveDrawWindow(sprint("Graphs_v2/Forecasts_CBS_scenario", j+1, ".pdf"));
			CloseDrawWindow();
			s_igr=0;
		}*/ 
		
		/*for (j=0; j<=3; j++)
		{
			for (i=0; i<=s_cMod-1; i++)
			{
				DrawTMatrix(s_igr, myforse[j*s_cMod+i][],  sprint("Forecast s.e. model ",  i+1), 2011, s_iFirstQuarter, 4, 0, i+2);
				DrawAdjust(ADJ_COLOR, i+2, 1);
			}
			SaveDrawWindow(sprint("Graphs_v2/Forecasts_uncertainty_CBS_scenario", j+1, ".pdf"));
			CloseDrawWindow();
		}*/ 

		// combined forecast and uncertainty graph for model 3, scenario 3
		/*s_igr=0;
		DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, s_iFirstQuarter, 4, 0, 1);
		DrawTMatrix(s_igr, myfor[2*s_cMod+1][], "Forecast model 2", 2011, s_iFirstQuarter, 4, 0, 2);
		DrawZ(myforse[2*s_cMod+1][], "95% CI", ZMODE_BAND, 1.96);
		DrawAdjust(ADJ_MINMAX, 61, 73);
		s_igr++;
		DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, s_iFirstQuarter, 4, 0, 1);
		DrawTMatrix(s_igr, myforCovid[2*s_cMod+2][], "Forecast model 3", 2011, s_iFirstQuarter, 4, 0, 2);
		DrawZ(myforseCovid[2*s_cMod+2][], "95% CI", ZMODE_BAND, 1.96);
		DrawAdjust(ADJ_MINMAX, 61, 73);
		DrawAdjust(ADJ_AREAMATRIX, 2, 1);
		SaveDrawWindow("Graphs_v2/Forecasts_CBS_winners.pdf");
		CloseDrawWindow();*/

	
}
