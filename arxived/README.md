# boundedcoal
In the folder boundsyndata (10 sequences), we generate synthetic datasets from bounded coalescent models.
syn1 works well, however, syn2 and syn2_2 doesn't work well. Not sure the reason currently.


In the folder standardsyndata (a single sequence), we mainly test for the case of syn2. We generate datasets based on standard coalscent model of syn2's Ne(t) and debug the code. Everything is under standard coalescent models.

Here's the summary for Ne(t):

syn1: Ne(t)=1

syn2: Ne(t)=25exp(-5t)

syn2_2: Ne(t)=25exp(-5t)+1
