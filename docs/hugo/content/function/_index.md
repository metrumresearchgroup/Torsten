+++
title = "Using Torsten"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-28T19:35:14-07:00
draft = false
weight = 1005
+++

<a id="org125f9fd"></a>

The reader should have a basic understanding of how Stan works before
reading this chapter. There are excellent resources online to get
started with Stan ([http://mc-stan.org/documentation](http://mc-stan.org/documentation)).
In this section we go through the different functions Torsten adds to
Stan. The code for the examples can be found at the `example-models` folder.

Torsten's functions are prefixed with `pmx_`.
For some of their arguments we adopt NM-TRAN format for events
specification(Table [tab:event_args](#tab:event_args)).

<a id="table--tab:event-args"></a>
<div class="table-caption">
  <span class="table-number"><a href="#table--tab:event-args">Table 1</a></span>:
  NM-TRAN compatible event specification arguments. All arrays should have the same length corresponding to the number of events.
</div>

| Argument Name | Definition                  | Stan data type |
|---------------|-----------------------------|----------------|
| `time`        | event time                  | `real[]`       |
| `amt`         | dosing amount               | `real[]`       |
| `rate`        | infusion rate               | `real[]`       |
| `ii`          | interdose interval          | `real[]`       |
| `evid`        | event ID                    | `int[]`        |
| `cmt`         | event compartment           | `int[]`        |
| `addl`        | additionial identical doses | `int[]`        |
| `ss`          | steady-state dosing flag    | `int[]`        |

All the `real[]` arguments above are allowed to
be `parameters` in a Stan model.
In addtion, Torsten functions
support optional arguments and overloaded signatures.
Optional arguments are indicated by surrounding square bracket `[]`.
Table below shows three commonly used PMX model arguments that support
overloading. In the rest of this document we assume this convention unless indicated otherwise.

<a id="table--tab:event-params"></a>
<div class="table-caption">
  <span class="table-number"><a href="#table--tab:event-params">Table 2</a></span>:
  PMX model parameter overloadings. One can use 1d array <code class="src src-stan"><span style="color: #b58900;">real</span>[]</code> to indicate constants of all events, or 2d array <code class="src src-stan"><span style="color: #b58900;">real</span>[ , ]</code> so that the \(i\)th row of the array describes the model arguments for time interval \((t_{i-1}, t_i)\), and the number of the rows equals to the size of <code>time</code>.
</div>

| Argument Name | Definition               | Stan data type          | Optional           |
|---------------|--------------------------|-------------------------|--------------------|
| `theta`       | model parameters         | `real[]` or `real[ , ]` | N                  |
| `biovar`      | bioavailability fraction | `real[]` or `real[ , ]` | Y (default to 1.0) |
| `tlag`        | lag time                 | `real[]` or `real[ , ]` | Y (default to 0.0) |
