Calibrating 2PL model via (a,b), proposed by Birnbaum (1968). Adjustment with Marginal Bayesian Modal Item Parameter Estimation.
Description
This function can estimate the item parameters of the 2PL model via Bayesian statistical methodology include Cornfield (1969), 
de Finetti (1974), Edwards, Lindman, and Savage (1963), Lindley (1970a, 1970b, 1971), and Novick and Jackson (1974). 
Examples of the use of these methods in educational settings can be found in Novick and Jackson (1974), Novick, Jackson, Thayer, 
and Cole (1972), and Rubin (1980). Lord (1986) compared maximum likelihood and Bayesian estimation methods in IRT.

Details
Two parameter logistic (2PL) model proposed by Birnbaum (1968):

P(x=1∣θ,a,b)=(1+exp(−D∗a∗(θ−b)))

where x=1 is the correct response, theta is examinne's ability; a, b are the item discrimination, and difficulty parameter, respectively;
D is the scaling constant 1 normal.

The implementation of expected a posteriori (EAP) ability estimation has communalities with the E-step of 
both marginal maximum likelihood / expectation-maximization (MMLE/EM) and Bayesian modal estimation / expectation-maximization 
(MBE/EM) estimation of item parameters.
