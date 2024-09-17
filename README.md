# boundedcoal
In the folder boundsyndata, we generate synthetic datasets from bounded coalescent models.
syn1 and syn2_2 works well, however, syn2 doesn't work well due to the right tails of Ne(t) approximating zero.


In the folder standardsyndata, we mainly test for the case of syn2. We generate datasets based on standard coalscent model of syn2's Ne(t) and debug the code. Everything is under standard coalescent models.

Here's the summary for Ne(t):
syn1: Ne(t)=1

syn2: Ne(t)=25exp(-5t)

syn2_2: Ne(t)=25exp(-5t)+1
