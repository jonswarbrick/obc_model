clear;

%% Fetch FRED data
clear;
c = fred('https://research.stlouisfed.org/fred2/');
freddata.gdp = fetch(c,'GDP'); % nominal gdp 1947-2015
freddata.rgdp = fetch(c,'GDPC1'); % real gdp 1947-2015
freddata.inv = fetch(c,'FPI'); % Fixed Private Investment 1947-2015
freddata.R = fetch(c,'FEDFUNDS'); % effective fed funds rate 1954-2016
freddata.P = fetch(c,'CPIAUCSL'); % CPI all urban consumers, all items 1947-2015
freddata.P_gdp = fetch(c,'GDPDEF'); % GDP deflator 1947-2015
freddata.PI = fetch(c,'GPDICTPI'); % CPI all urban consumers, all items 1947-2015
freddata.spread = fetch(c,'BAA10YM'); % Moody's Seasoned Baa Corporate Bond Yield Relative to Yield on 10-Year Treasury Constant Maturity 1953-2016
freddata.AAAyield = fetch(c,'AAA'); % Moody's Seasoned Aaa Corporate Bond Yield 1953-2016
freddata.BAAyield = fetch(c,'BAA'); % Moody's Seasoned Baa Corporate Bond Yield 1953-2016
freddata.BAAFFM = fetch(c,'BAAFFM'); % Moody's Seasoned Baa Corporate Bond Yield to FFR spread 1953-2016
freddata.C = fetch(c,'PCEC'); % Personal Consumption Expenditures 1947-2015
freddata.av_H = fetch(c,'PRS85006023'); % Nonfarm Business Sector: Average Weekly Hours 1947-2015
freddata.Emp = fetch(c,'CE16OV'); % Civilian Employment 1948-2016
freddata.G = fetch(c,'GCE'); % Government Consumption Expenditures and Gross Investment  1947-2015
freddata.pop = fetch(c,'CNP16OV'); % Civilian Noninstitutional Population  1948-2016

save('rawData.mat','freddata')
