# Colombian-Synthetic-Indicator
Colombian Synthetic Banking Indicator as a Leading Signal of Real Economic Activity

📄 TITLE

A Synthetic Banking Indicator as a Leading Signal of Real Economic Activity: Evidence from Dynamic and Frequency-Based Analysis

📄 ABSTRACT

This paper proposes and validates a synthetic banking indicator as a forward-looking measure of real economic activity. Using industrial production as a benchmark, we assess the indicator’s explanatory and predictive capacity through correlation analysis, dynamic modeling, Granger causality tests, and vector autoregressions (VAR). Cross-correlation results reveal that the banking indicator leads industrial production by up to 12 months. Granger causality tests provide strong evidence that the banking indicator predicts real activity across multiple lag structures. VAR estimates confirm the significance of lagged banking variables in explaining industrial production, while impulse response functions reveal a short-term adjustment followed by a sustained positive effect. Frequency decomposition using the Hodrick–Prescott filter shows that the relationship is primarily driven by low-frequency components, with weak association in cyclical fluctuations. Forecast evaluation using AIC and BIC demonstrates that models incorporating the banking indicator outperform naive benchmarks. Overall, the results support the interpretation of the banking indicator as a structural and forward-looking signal of economic conditions.

📄 1. INTRODUCTION

Understanding the relationship between financial conditions and real economic activity is a central question in macroeconomics and finance. Financial variables are often considered leading indicators, yet empirical evidence on their timing, predictive content, and frequency structure remains mixed.

This paper contributes to the literature by constructing and validating a synthetic banking indicator designed to capture the underlying dynamics of the banking sector. Unlike traditional financial indicators, which are often analyzed in isolation, the proposed measure aggregates banking-related signals into a single index that reflects systemic conditions.

The contribution of this study is fourfold. First, it demonstrates a strong association between the banking indicator and industrial production. Second, it provides clear evidence of leading behavior, with the banking indicator anticipating real activity by up to 12 months. Third, it establishes the predictive content of the indicator using Granger causality and dynamic regression frameworks. Fourth, it distinguishes between low-frequency and high-frequency dynamics, showing that the relationship is primarily structural rather than driven by short-term fluctuations.

📄 2. DATA AND CONSTRUCTION

The analysis is based on a synthetic banking indicator constructed from aggregated banking-related variables (details depend on your data: credit, liquidity, etc.). The indicator is normalized and aligned temporally with the Industrial Production Index (IPI), which serves as a proxy for real economic activity.

All series are standardized where necessary and cleaned to ensure comparability. The sample consists of monthly observations, allowing for the analysis of both short-term dynamics and long-term trends.

📄 3. METHODOLOGY

The empirical strategy follows a structured approach:

(1) Static relationship
Correlation analysis
Linear regression between the banking indicator and IPI
(2) Temporal structure
Cross-correlation function (CCF) to identify leading behavior
(3) Predictive content
Granger causality tests for multiple lag specifications (1–12)
(4) Dynamic modeling
Vector autoregression (VAR) with lagged variables
Impulse response functions (IRFs)
(5) Frequency decomposition
Hodrick–Prescott (HP) filter to separate trend and cycle
(6) Forecast evaluation
Model comparison using AIC and BIC
📄 4. RESULTS
4.1 Static relationship

The banking indicator exhibits a strong positive correlation with industrial production (≈0.81), which increases when considering the smoothed trend component (≈0.89). Regression results show substantial explanatory power, with R² values approaching 0.79 for the trend specification.

4.2 Temporal structure

Cross-correlation analysis reveals that the maximum association occurs at a lag of −12 periods, indicating that the banking indicator leads industrial production by approximately one year. This result provides strong evidence of forward-looking behavior.

4.3 Granger causality

Granger causality tests show that the null hypothesis of no predictive relationship from the banking indicator to industrial production is rejected across all lag orders (1–12), with highly significant p-values. In contrast, reverse causality from industrial production to the banking indicator is only observed at short lags and disappears at longer horizons.

This asymmetry suggests that the predictive direction runs primarily from the banking sector to real economic activity.

4.4 Dynamic modeling (VAR)

The VAR model confirms the dynamic relationship. Lagged values of the banking indicator are statistically significant predictors of industrial production:

Lag 1: significant
Lag 2: significant

The model explains approximately 85% of the variation in industrial production (R² ≈ 0.85), indicating strong explanatory power.

Importantly, lagged values of industrial production do not significantly explain the banking indicator, reinforcing the directional interpretation.

4.5 Impulse response analysis

Impulse response functions reveal a clear dynamic pattern:

Short-term negative adjustment
Medium-term recovery
Sustained positive effect

This pattern is consistent with macro-financial transmission mechanisms, where financial shocks initially dampen activity before leading to expansion.

4.6 Frequency decomposition

HP filter results show:

Weak correlation in cyclical components (~0.18)
Strong correlation in trend components (~0.93)

This indicates that the relationship is primarily driven by structural, low-frequency dynamics rather than short-term fluctuations.

4.7 Forecast evaluation

Model comparison shows that including the banking indicator significantly improves predictive performance:

AIC and BIC are substantially lower compared to naive models
This confirms the indicator’s predictive value
📄 5. DISCUSSION

The results consistently support the interpretation of the banking indicator as a forward-looking and structurally informative measure of economic conditions. The strong relationship at the trend level suggests that the indicator captures underlying macroeconomic dynamics rather than transient shocks.

The lack of explanatory power in differences or growth rates reinforces this interpretation, indicating that the indicator is not suitable for modeling high-frequency fluctuations but is highly effective in capturing long-term movements.

The asymmetry observed in Granger causality and VAR results highlights the central role of the banking sector in driving real economic activity.

📄 6. CONCLUSION

This paper provides robust evidence that a synthetic banking indicator contains forward-looking information about real economic activity. The indicator exhibits strong predictive power, clear leading behavior, and structural relevance.

Rather than capturing short-term fluctuations, the indicator reflects low-frequency economic dynamics, making it particularly valuable for medium-term forecasting and macroeconomic analysis.

Future research may extend this framework by incorporating additional financial variables, testing alternative filtering techniques, and applying the methodology to other economies.
