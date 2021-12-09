# Discussion on prior elicitation

# Prior elicitation can be as simple as drawing a PDF
curve(dnorm(x, mean=0, sd=1), from=-4, to=4)
curve(dnorm(x, mean=3.5, sd=1), from=0, to=8)
curve(dnorm(x, mean=3.5, sd=4), from=-10, to=20)

# e.g. probability parameter
curve(dunif(x, 0, 1))
curve(dbeta(x, 2, 1.2))
