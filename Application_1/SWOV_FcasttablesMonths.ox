#include <oxstd.oxh>
#include <oxdraw.h>

printlatextable(const file, const stitle, const asformat, const m)
{
	decl i, j, asmonths, ncol = columns(m), nrow = rows(m);
	decl s;
	s = sprint("\\","begin{table}[ht]\n");
	s = sprint(s, "\\","centering\n");
	s = sprint(s, "\\","caption{", stitle, "}\n");
	s = sprint(s, "\\","resizebox{\columnwidth}{!}{\n");
	s = sprint(s, "\\","begin{tabular}{l");
	for (i=0; i<ncol; i++) s = sprint(s, "c");
	s = sprint(s, "}\n");
	s = sprint(s, "\\","hline\n");
	asmonths = {"Jan", "Feb", "Mar", "Apr", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec"};
	s = sprint(s, " & ", "Model");
	for (i=0; i<sizec(asmonths); i++) s=sprint(s, " &", "%5d", asmonths[i]);
	//s = sprint(s, " & ", "Total");
	s = sprint(s, " ", "\\", "\\", "\n");
	s = sprint(s, "\\","hline\n");			  
	for (j=0; j<nrow; j++)
	{
		for (i=0; i<ncol; i++)
			s = (sizeof(asformat)>0) ? sprint(s, " & ", asformat[i], m[j][i]) : sprint(s, " & ", m[j][i]);
		s = sprint(s, " ", "\\", "\\", "\n");
	}
	s = sprint(s, "\\","hline\n");
	s = sprint(s, "\\","end{tabular}}\n");
	s = sprint(s, "\\","end{table}\n");

	println(s);
	if (isfile(file))
		fprintln(file, s);
}


main()
{
	decl mdata = log(loadmat("Gecorrigeerd_Aangevuld_NewMaandcijfers_swov.csv"))';
	println("mdata", mdata);
	decl vyobstot = mdata[1][290:433];
	decl myfor = loadmat("Output_v2/myfor.mat"), myforse = loadmat("Output_v2/myforse.mat");
	decl myforPRERFDC = loadmat("Output_v2/myforPRERFDC.mat"), myforsePRERFDC = loadmat("Output_v2/myforsePRERFDC.mat");
	decl myforFINRFPRERFDC = loadmat("Output_v3/myforFINRF-PRERFDC.mat"), myforseFINRFPRERFDC = loadmat("Output_v3/myforseFINRF-PRERFDC.mat");
	decl s_cMod = 6, s_igr = 0;
	decl s_sNameY = {"FINRF", "REGRF", "PRERFDC", "FINRFDC"};
	println("vyobstot size: ", columns(vyobstot));
	println("vyobstot: ", vyobstot);

	decl i, j=3, merr, merrPRERFDC, merrFINRFPRERFDC;
	merr = vyobstot-myfor;
	println(merr);
	merrPRERFDC = vyobstot-myforPRERFDC;
	println("merrPRERFDC", merrPRERFDC);
	merrFINRFPRERFDC = vyobstot-myforFINRFPRERFDC;
	println("merrFINRFPRERFDC", merrFINRFPRERFDC);

	decl mMSE = <>, mMAE = <>, mMAPE = <>, mMSEPRERFDC = <>, mMSEFINRFPRERFDC = <>;
	for (i=0; i<12; i++)
	{
		decl mMonth = <>, mMonthrel = <>, mMonthPRERFDC = <>, mMonthFINRFPRERFDC = <>;
		for (j=0; j<144; j=j+12)
		{
			mMonth ~= merr[][j+i];
			mMonthPRERFDC ~= merrPRERFDC[][j+i];
			mMonthFINRFPRERFDC ~= merrFINRFPRERFDC[][j+i];
			mMonthrel ~= merr[][j+i]./vyobstot[j+i];
		}
		println("mMonth", mMonth);
		mMSE ~= meanr(sqr(mMonth));
		mMSEPRERFDC ~= meanr(sqr(mMonthPRERFDC));
		mMSEFINRFPRERFDC ~= meanr(sqr(mMonthFINRFPRERFDC));
		//println("mMSE", mMSE);
		mMAE ~= meanr(fabs(mMonth));
		//println("mMAE", mMAE);
		mMAPE ~= meanr(100*fabs(mMonthrel));
	}
	println("MSE", mMSE);

	decl vMSEbase=mMSE[0][], vMAEbase=mMAE[0][], vMAPEbase=mMAPE[0][];
	mMSE ./= vMSEbase;
	mMAE ./= vMAEbase;
	mMAPE ./= vMAPEbase;
	mMSEPRERFDC ./= vMSEbase;
	mMSEFINRFPRERFDC ./= vMSEbase;

	decl file_mse_months = fopen("table_mse_months.tex", "w"), asformat = new array[14];
	asformat[0] = sprint("%5d");
	for (i=1; i<=13; i++) asformat[i] = sprint("%7.3f");
	printlatextable(file_mse_months, sprint("MSE forecasting performance based on the 1- to 12-step ahead forecasts per month."), asformat, ones(5,1) ** range(1, 6)' ~ mMSE);
	printlatextable(<>, sprint("MSE forecasting performance based on the 1- to 12-step ahead forecasts per year and in total of univariate stsm of PRERFDC."), asformat, ones(4,1) ** range(1, 2)' ~ mMSEPRERFDC);
	printlatextable(<>, sprint("MSE forecasting performance based on the 1- to 12-step ahead forecasts per year and in total of bivariate stsm of FINRF and PRERFDC."), asformat, ones(4,1) ** range(1, 2)' ~ mMSEFINRFPRERFDC);
	printlatextable(<>, sprint("MAE forecasting performance based on the 1- to 12-step ahead forecasts per month."), asformat, ones(5,1) ** range(1, 6)' ~ mMAE);
	printlatextable(<>, sprint("MAPE forecasting performance based on the 1- to 12-step ahead forecasts per month."), asformat, ones(5,1) ** range(1, 6)' ~ mMAPE);

}
