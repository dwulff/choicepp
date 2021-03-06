% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ITCH_lik}
\alias{ITCH_lik}
\title{ITCH deterministic choices}
\usage{
ITCH_lik(par, problems, choices, epsilon = 1e-04, type = 0L)
}
\arguments{
\item{par}{numeric vector specifying the discounting factor
and the choice sensitivity (in that order).}

\item{problems}{a numeric matrix with four columns containing for each decision
the outcome and delay of the sooner smaller option, the outcome and delay
of the larger later option (in that order).}

\item{choices}{numeric vector containing the choice choice (0 = sooner smaller,
1 = larger later).}

\item{epsilon}{numeric scalar specifying a boundary on the likelihood of a choice.
Specifically, likelihoods will be bound between epsilon and 1-epsilon.}

\item{type}{integer specifying the model variant to be used. See details.}
}
\value{
The negative log-liklihood.
}
\description{
Generate deterministic choice from various intertemporal choice heuristics (ITCH).
}
\details{
Function takes parameter values and the choice problem to generate
deterministic choices.

Functions implements different variants as a function of \code{type}:

00 implements the four parameter variant with exponential choice rule

10 implements the five parameter variant with exponential choice rule

The four and five parameter variants are defined as
\deqn{\beta_0 * (o_ss - o_ll) + \beta_1 * \frac{o_ss - o_ll}{o*} +
      \beta_2 * (d_ss - d_ll) + \beta_3 * \frac{d_ss - d_ll}{d*}}
and

\deqn{\beta_5 +
      \beta_0 * (o_ss - o_ll) + \beta_1 * \frac{o_ss - o_ll}{o*} +
      \beta_2 * (d_ss - d_ll) + \beta_3 * \frac{d_ss - d_ll}{d*}}

With o_ss, o_ll, d_ss, d_ll being the outcomes and delays of the
sooner smaller and larger later options, respectively, and o* and d*
are the arithmetic means of both outcomes and delays, respectively.

Note that the ITCH cannot be used with a homothetic choice rule.
}
