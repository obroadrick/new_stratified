"""
Script to mess around and hopefully
develop some intuition
for convolutions.
"""

# ROUND 1
n1 = 100
p = .7
# normal binomial dist
dist = binom.pmf(range(0,n1+1),n1,p)

# DRAW FOR ROUND 1
k = 80
# chop off tail
for i in range(80,n+1):
    dist[i] = 0

# ROUND 2
n2 = 60
# new binomial dist
dist = binom.pmf(range(0,n2+1),n2,p)
# convolve
for i in range(0,n1+1):
    for j in range(0,n2+1):
        dist[i]


"""
AH OK SO AT THIS POINT IN WRITING
THE SCRIPT I REALIZED WHAT I WAS
CONFUSED ABOUT AND HAVE MY UNDERSTANDING
SORTED OUT AGAIN... CONVOLUTIONS
ARE NEEDED SINCE SOME OF THE WAYS
OF OBTAINING LESS THAN THE IMPOSSIBLE
NUMBERS IN THE CONVOLVED DISTRIBUTION
ARE REMOVED NOT JUST THE ONES THAT ARE
NOT POSSIBLE ANY MORE
"""


