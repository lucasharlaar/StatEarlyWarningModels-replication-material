#include <oxstd.oxh>
#include <oxdraw.h>

main()
{
	decl mdata = log(loadmat("Gecorrigeerd_Aangevuld_NewMaandcijfers_swov.csv"))';
	decl vyobstot = mdata[1][290:433];
	decl myfor = loadmat("Output_v2/myfor.mat"), myforse = loadmat("Output_v2/myforse.mat");
	decl s_cMod = 6, s_igr = 0;
	decl s_sNameY = {"FINRF", "REGRF", "PRERFDC", "FINRFDC"};
	println("vyobstot size: ", columns(vyobstot));
	println("vyobstot: ", exp(vyobstot));

	decl i, j=3;
		// scenario comparison plots
		/*for (j=0; j<=s_cMod-1; j++)
		{
			s_igr=0;
			DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, 1, 12, 0, 1);
			DrawAdjust(ADJ_COLOR, 1, 14);
			DrawTMatrix(s_igr, myfor[2][], "Forecast univariate", 2011, 1, 12, 0);
			DrawAdjust(ADJ_COLOR, 2, 9);
			for (i=1; i<=2; i++)
			{
				DrawTMatrix(s_igr, myfor[i*s_cMod+j][],  sprint("Forecast in scenario ",  i+1), 2011, 1, 12, 0);
				DrawAdjust(ADJ_COLOR, i+2, 7);
				DrawAdjust(ADJ_MINMAX, 3.2, 4.7);
				//DrawText(s_igr, "log(RF)", 0, 0, -1, 300, TEXT_YLABEL);
				//DrawText(s_igr, "Time (month)", 0, 0, -1, 300, TEXT_XLABEL);
			}
			s_igr++;
			DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, 1, 12, 0, 1);
			DrawAdjust(ADJ_COLOR, 1, 14);
			DrawTMatrix(s_igr, myfor[2][], "Forecast univariate", 2011, 1, 12, 0);
			DrawAdjust(ADJ_COLOR, 2, 9);
			for (i=3; i<=4; i++)
			{
				DrawTMatrix(s_igr, myfor[i*s_cMod+j][],  sprint("Forecast in scenario ",  i+1), 2011, 1, 12, 0);
				DrawAdjust(ADJ_COLOR, i+2, 7);
				DrawAdjust(ADJ_MINMAX, 3.2, 4.7);
				//DrawText(s_igr, "log(RF)", 0, 0, -1, 300, TEXT_YLABEL);
				//DrawText(s_igr, "Time (month)", 0, 0, -1, 300, TEXT_XLABEL);
			}
			DrawAdjust(ADJ_AREAMATRIX, 2,1);
			SaveDrawWindow(sprint("Graphs_v2/ForecastsSWOV_model", j+1, ".pdf"));
			CloseDrawWindow();
		}*/
		/*
		//for (j=0; j<=s_cMod-1; j++)
		//{
			DrawTMatrix(s_igr, myforse[j][], "Forecast s.e. univariate", 2011, 1, 12, 0);
			DrawAdjust(ADJ_COLOR, 2, 9);
			for (i=1; i<=4; i++)
			{
				DrawTMatrix(s_igr, myforse[i*s_cMod+j][],  sprint("Forecast s.e. in scenario ",  i+1), 2011, 1, 12, 0);
				DrawAdjust(ADJ_COLOR, i+2, 7);
				DrawAdjust(ADJ_MINMAX, 0.025, 0.25);
			}
			//SaveDrawWindow(sprint("Graphs_v2/Forecasts_uncertaintySWOV_model", j+1, ".pdf"));
			//CloseDrawWindow();
		//}*/
		// model comparison plots
		for (j=0; j<=4; j++)
		{
			s_igr=0;
			DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, 1, 12, 0, 1);
			DrawAdjust(ADJ_COLOR, 1, 14);
			DrawTMatrix(s_igr, myfor[0][], "Forecast univariate", 2011, 1, 12, 0);
			DrawAdjust(ADJ_COLOR, 2, 9);
			for (i=0; i<=2; i++)
			{
				DrawTMatrix(s_igr, myfor[j*s_cMod+i][],  sprint("Forecast model ",  i+1), 2011, 1, 12, 0);
				DrawAdjust(ADJ_COLOR, i+3, 7);
				DrawAdjust(ADJ_MINMAX, 3.2, 4.7);
			}
			s_igr++;
			DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, 1, 12, 0, 1);
			DrawAdjust(ADJ_COLOR, 1, 14);
			DrawTMatrix(s_igr, myfor[1][], "Forecast univariate", 2011, 1, 12, 0);
			DrawAdjust(ADJ_COLOR, 2, 9);
			for (i=3; i<=s_cMod-1; i++)
			{
				DrawTMatrix(s_igr, myfor[j*s_cMod+i][],  sprint("Forecast model ",  i+1), 2011, 1, 12, 0);
				DrawAdjust(ADJ_COLOR, i+3, 7);
				DrawAdjust(ADJ_MINMAX, 3.2, 4.7);
			}
			DrawAdjust(ADJ_AREAMATRIX, 2, 1);
			SaveDrawWindow(sprint("Graphs_v2/ForecastsSWOV_scenario", j+1, ".pdf"));
			CloseDrawWindow();
		}
		/*
		//for (j=0; j<=4; j++)
		//{
		j=2;
		s_igr++;
			DrawTMatrix(s_igr, myforse[0][], "Forecast s.e. univariate", 2011, 1, 12, 0);
			DrawAdjust(ADJ_COLOR, 2, 9);
			for (i=0; i<=s_cMod-1; i++)
			{
				DrawTMatrix(s_igr, myforse[j*s_cMod+i][],  sprint("Forecast s.e. model ",  i+1), 2011, 1, 12, 0);
				DrawAdjust(ADJ_COLOR, i+3, 7);
				DrawAdjust(ADJ_MINMAX, 0.025, 0.25);
			}
			//SaveDrawWindow(sprint("Graphs_v2/Forecasts_uncertaintySWOV_scenario", j+1, ".pdf"));
			//CloseDrawWindow();
		//}
			DrawAdjust(ADJ_AREAMATRIX, 2, 1);
			//DrawText(s_igr+1, "log(RF)", 0, 0, -1, 250, TEXT_YLABEL);
			//DrawText(s_igr+1, "Time (month)", 0, 0, -1, 250, TEXT_XLABEL);
			//DrawAxisAuto(s_igr+1,0,FALSE);
			//DrawAxisAuto(s_igr+1,1,FALSE);
			SaveDrawWindow(sprint("Graphs_v2/Forecasts_uncertaintySWOV_mod4scen3.pdf"));
			CloseDrawWindow();
		//}
		*/
		// combined forecast and uncertainty graphs for winning models in scenario 3
		/*DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, 1, 12, 0, 1);
		DrawTMatrix(s_igr, myfor[2*s_cMod+4][], "Forecast model 5", 2011, 1, 12, 0, 2);
		decl vCI_L, vCI_U;
		vCI_L = myfor[2*s_cMod+4][]-1.96*myforse[2*s_cMod+4][];
		vCI_U = myfor[2*s_cMod+4][]+1.96*myforse[2*s_cMod+4][];
		DrawTMatrix(s_igr, vCI_L, "95% CI", 2011, 1, 12, 0);
		DrawAdjust(ADJ_COLOR, 2, 5);
		DrawTMatrix(s_igr, vCI_U, "", 2011, 1, 12, 0);
		DrawAdjust(ADJ_COLOR, 2, 5);
		//DrawZ(myforse[2*s_cMod+4][], "95% CI", ZMODE_BAND, 1.96);
		//DrawAdjust(ADJ_COLOR, 8, 3);
		DrawAdjust(ADJ_MINMAX, 3.2, 4.7);
		s_igr++;
		DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, 1, 12, 0, 1);
		DrawTMatrix(s_igr, myfor[4*s_cMod+4][], "Forecast model 5", 2011, 1, 12, 0, 2);
		vCI_L = myfor[4*s_cMod+4][]-1.96*myforse[4*s_cMod+4][];
		vCI_U = myfor[4*s_cMod+4][]+1.96*myforse[4*s_cMod+4][];
		DrawTMatrix(s_igr, vCI_L, "95% CI", 2011, 1, 12, 0);
		DrawAdjust(ADJ_COLOR, 2, 5);
		DrawTMatrix(s_igr, vCI_U, "", 2011, 1, 12, 0);
		DrawAdjust(ADJ_COLOR, 2, 5);
		DrawAdjust(ADJ_MINMAX, 3.2, 4.7);
		DrawAdjust(ADJ_AREAMATRIX, 2, 1);
		SaveDrawWindow("Graphs_v2/Forecast_model5_scen3_log.pdf");
		CloseDrawWindow();*/
		/*DrawTMatrix(s_igr, vyobstot,  s_sNameY[0], 2011, 1, 12, 0);
		DrawAdjust(ADJ_COLOR, 1, 2);
		DrawTMatrix(s_igr, myfor[2*s_cMod+5][], "Forecast in scenario 3", 2011, 1, 12, 0, 2);
		DrawZ(myforse[2*s_cMod+5][], "", ZMODE_BAND, 1.96);
		DrawAdjust(ADJ_COLOR, 2, 5);
		SaveDrawWindow("Graphs_v2/Forecast_model6_scen3_log.pdf");
		CloseDrawWindow();*/
		/*DrawTMatrix(s_igr, exp(vyobstot),  s_sNameY[0], 2011, 1, 12, 0, 1);
		//DrawAdjust(ADJ_COLOR, 1, 2);
		DrawTMatrix(s_igr, exp(myfor[2*s_cMod+4][]), "Forecast model 5", 2011, 1, 12, 0, 2);
		//DrawAdjust(ADJ_COLOR, 8, 3);
		SaveDrawWindow("Graphs_v2/Forecast_model5_scen3.pdf");
		CloseDrawWindow();
		/*DrawTMatrix(s_igr, exp(vyobstot),  s_sNameY[0], 2011, 1, 12, 0);
		DrawAdjust(ADJ_COLOR, 1, 2);
		DrawTMatrix(s_igr, exp(myfor[2*s_cMod+5][]), "Forecast in scenario 3", 2011, 1, 12, 0, 2);
		SaveDrawWindow("Graphs_v2/Forecast_model6_scen3.pdf");
		CloseDrawWindow();*/  */
	
}
