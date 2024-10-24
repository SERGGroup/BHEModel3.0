# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from pymc import HalfCauchy, Model, Normal, sample
import matplotlib.pyplot as plt
import pandas as pd
import arviz as az
import numpy as np
import pymc as pm

print(f"Running on PyMC v{pm.__version__}")

# %%------------   GENERATE RANDOM DATA                   -----------------------------------------------------------> #
RANDOM_SEED = 8927
rng = np.random.default_rng(RANDOM_SEED)

size = 200
true_intercept = 1
true_slope = 2

x = np.linspace(0, 1, size)
# y = a + b*x
true_regression_line = true_intercept + true_slope * x
# add noise
y = true_regression_line + rng.normal(scale=0.5, size=size)

data = pd.DataFrame(dict(x=x, y=y))

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, xlabel="x", ylabel="y", title="Generated data and underlying model")
ax.plot(x, y, "x", label="sampled data")
ax.plot(x, true_regression_line, label="true regression line", lw=2.0)
plt.legend(loc=0)
plt.show()


# %%------------   ESTIMATING THE MODEL                   -----------------------------------------------------------> #
with Model() as model:  # model specifications in PyMC are wrapped in a with-statement
    # Define priors
    sigma = HalfCauchy("sigma", beta=10)
    intercept = Normal("Intercept", 0, sigma=20)
    slope = Normal("slope", 0, sigma=20)

    # Define likelihood
    likelihood = Normal("y", mu=intercept + slope * x, sigma=sigma, observed=y)

    # Inference!
    # draw 3000 posterior samples using NUTS sampling
    idata = sample(3000)


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
az.plot_trace(idata, figsize=(10, 7))
plt.tight_layout()
plt.show()
