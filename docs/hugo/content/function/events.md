+++
title = "Events specification"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-30T11:38:23-07:00
draft = false
weight = 2001
+++

Torsten's functions are prefixed with `pmx_`.
For some of their arguments we adopt NM-TRAN format for events
specification(Table [1](#table--tab:event-args)).

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
