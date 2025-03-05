#!/bin/sh

set -eu

# Check that stanc3 recognizes each of the entry points below by
# injecting the entry point into a dummy model and verifying that
# stanc3 doesn't complain about an undeclared identifier.

entrypoints='
pmx_ln_interpolate
pmx_ode_adams
pmx_ode_adams_ctrl
pmx_ode_bdf
pmx_ode_bdf_ctrl
pmx_ode_ckrk
pmx_ode_ckrk_ctrl
pmx_ode_rk45
pmx_ode_rk45_ctrl
pmx_solve_adams
pmx_solve_bdf
pmx_solve_group_adams
pmx_solve_group_bdf
pmx_solve_group_rk45
pmx_solve_linode
pmx_solve_onecpt
pmx_solve_onecpt_bdf
pmx_solve_onecpt_effcpt
pmx_solve_onecpt_rk45
pmx_solve_rk45
pmx_solve_twocpt
pmx_solve_twocpt_bdf
pmx_solve_twocpt_effcpt
pmx_solve_twocpt_rk45
'

test -d stanc3 || {
    printf >&2 '%s must be called from top of Torsten repository\n'  "$0"
    exit 2
}

STANC3=$PWD/stanc3
export STANC3

command -v ocaml >/dev/null || {
    printf >&2 'ocaml not in PATH; is the stanc3 dev environment activated?\n'
    exit 2
}

tdir=$(mktemp -d "${TMPDIR:-/tmp}"/torsten-check-entrypoints-XXXXX)
trap 'rm -rf "$tdir"' 0

printf >&2 'building cmdstan...\n'
make -C cmdstan -j "$(nproc)" build >/dev/null 2>"$tdir/build-stderr" || {
    cat >&2 "$tdir/build-stderr"
    exit 2
}

build_model () {
    cat<<EOF >"$tdir"/dummy.stan
transformed parameters {
 real x = $1();
}
EOF
    # shellcheck disable=SC2069
    make -C cmdstan "$tdir"/dummy 2>&1 >/dev/null
}

IFS='
'
failed=
for e in $entrypoints
do
    printf >&2 '%s...'  "$e"
    if build_model "$e" | grep -q 'undeclared identifier'
    then
        printf 'FAILED\n'
        failed=1
    else
        printf 'ok\n'
    fi
done
test -z "$failed"
