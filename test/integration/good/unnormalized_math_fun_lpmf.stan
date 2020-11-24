functions {
    real foo_lpmf(int i1, real r1, matrix m1, vector v1, int[] ai1, row_vector rv1) {
        real r;
        r += bernoulli_logit_glm_lupmf(i1| m1, r1, v1);
        r += bernoulli_logit_lupmf(i1| r1);
        r += bernoulli_lupmf(i1| r1);
        r += beta_binomial_lupmf(i1| i1, r1, r1);
        r += binomial_logit_lupmf(i1| i1, r1);
        r += binomial_lupmf(i1| i1, r1);
        r += categorical_logit_glm_lupmf(i1| rv1, v1, m1);
        r += categorical_logit_lupmf(i1| v1);
        r += categorical_lupmf(i1| v1);
        r += hypergeometric_lupmf(i1| i1, i1, i1);
        r += multinomial_lupmf(ai1| v1);
        r += neg_binomial_2_log_glm_lupmf(i1| m1, r1, v1, r1);
        r += neg_binomial_2_log_lupmf(i1| r1, r1);
        r += neg_binomial_2_lupmf(i1| r1, r1);
        r += neg_binomial_lupmf(i1| r1, r1);
        r += ordered_logistic_glm_lupmf(i1| rv1, v1, v1);
        r += ordered_logistic_lupmf(i1| r1, v1);
        r += ordered_probit_lupmf(i1| r1, v1);
        r += poisson_log_glm_lupmf(i1| m1, r1, v1);
        r += poisson_log_lupmf(i1| r1);
        r += poisson_lupmf(i1| r1);
        return r;
    }
    void goo_lp(int i1, real r1, matrix m1, vector v1, int[] ai1, row_vector rv1) {
        target += bernoulli_logit_glm_lupmf(i1| m1, r1, v1);
        target += bernoulli_logit_lupmf(i1| r1);
        target += bernoulli_lupmf(i1| r1);
        target += beta_binomial_lupmf(i1| i1, r1, r1);
        target += binomial_logit_lupmf(i1| i1, r1);
        target += binomial_lupmf(i1| i1, r1);
        target += categorical_logit_glm_lupmf(i1| rv1, v1, m1);
        target += categorical_logit_lupmf(i1| v1);
        target += categorical_lupmf(i1| v1);
        target += hypergeometric_lupmf(i1| i1, i1, i1);
        target += multinomial_lupmf(ai1| v1);
        target += neg_binomial_2_log_glm_lupmf(i1| m1, r1, v1, r1);
        target += neg_binomial_2_log_lupmf(i1| r1, r1);
        target += neg_binomial_2_lupmf(i1| r1, r1);
        target += neg_binomial_lupmf(i1| r1, r1);
        target += ordered_logistic_glm_lupmf(i1| rv1, v1, v1);
        target += ordered_logistic_lupmf(i1| r1, v1);
        target += ordered_probit_lupmf(i1| r1, v1);
        target += poisson_log_glm_lupmf(i1| m1, r1, v1);
        target += poisson_log_lupmf(i1| r1);
        target += poisson_lupmf(i1| r1);
    }
}
data {
    int i1;
    int ai1[5];
}
parameters {
    real r1;
    vector[5] v1;
    row_vector[5] rv1;
    matrix[5, 5] m1;
}
model {
    real r;
    r += bernoulli_logit_glm_lupmf(i1| m1, r1, v1);
    r += bernoulli_logit_lupmf(i1| r1);
    r += bernoulli_lupmf(i1| r1);
    r += beta_binomial_lupmf(i1| i1, r1, r1);
    r += binomial_logit_lupmf(i1| i1, r1);
    r += binomial_lupmf(i1| i1, r1);
    r += categorical_logit_glm_lupmf(i1| rv1, v1, m1);
    r += categorical_logit_lupmf(i1| v1);
    r += categorical_lupmf(i1| v1);
    r += hypergeometric_lupmf(i1| i1, i1, i1);
    r += multinomial_lupmf(ai1| v1);
    r += neg_binomial_2_log_glm_lupmf(i1| m1, r1, v1, r1);
    r += neg_binomial_2_log_lupmf(i1| r1, r1);
    r += neg_binomial_2_lupmf(i1| r1, r1);
    r += neg_binomial_lupmf(i1| r1, r1);
    r += ordered_logistic_glm_lupmf(i1| rv1, v1, v1);
    r += ordered_logistic_lupmf(i1| r1, v1);
    r += ordered_probit_lupmf(i1| r1, v1);
    r += poisson_log_glm_lupmf(i1| m1, r1, v1);
    r += poisson_log_lupmf(i1| r1);
    r += poisson_lupmf(i1| r1);
    r += foo_lupmf(i1| r1, m1, v1, ai1, rv1);
    goo_lp(i1, r1, m1, v1, ai1, rv1);
}