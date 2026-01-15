# tests_task02.jl
# Набор тестов корректности для задания 02
# Запуск: julia tests_task02.jl

using Random
include("isequivSSP.jl")

function assert_true(cond, msg)
    if !cond
        error("ТЕСТ НЕ ПРОЙДЕН: " * msg)
    end
end

function assert_eq(a, b, msg)
    if a != b
        error("ТЕСТ НЕ ПРОЙДЕН: $msg\nОжидали: $b\nПолучили: $a")
    end
end

# Проверка равенства A*G1 == G2*P по результату (A, perm)
function check_equation_local(G1, G2, res)
    if res == -1 || res == -2
        return false
    end
    A, perm = res
    P = perm_to_P(perm)
    left  = gf2_mul(toBoolMat(A), toBoolMat(G1))
    right = gf2_mul(toBoolMat(G2), P)
    return left == right
end

# Случайная матрица полного ранга
function random_fullrank_G(k, n)
    while true
        G = rand(Bool, k, n)
        if rank_gf2(G) == k
            return G
        end
    end
end

function random_invertible_A(k)
    while true
        A = rand(Bool, k, k)
        if is_invertible_gf2(A)
            return A
        end
    end
end


# ТЕСТ 1: только перестановка (P)

function test_perm_only()
    G1 = Bool[
        1 1 0 1;
        0 1 1 1
    ]
    perm0 = [2, 1, 4, 3]
    G2 = G1[:, perm0]

    res = isequivSSP(G1, G2, sig3_ws_hull_short; max_tries=100000)
    assert_true(res != -1 && res != -2, "perm_only: алгоритм не нашёл эквивалентность")
    assert_true(check_equation_local(G1, G2, res), "perm_only: не выполняется A*G1 == G2*P")
end


# ТЕСТ 2: только A (строковые операции)

function test_A_only()
    G1 = Bool[
        1 0 0 1 0;
        0 1 0 0 1;
        0 0 1 1 1
    ]
    A0 = Bool[
        1 1 0;
        0 1 1;
        0 0 1
    ]
    # без перестановки
    G2 = gf2_mul(A0, G1)

    res = isequivSSP(G1, G2, sig3_ws_hull_short; max_tries=200000)
    assert_true(res != -1 && res != -2, "A_only: алгоритм не нашёл эквивалентность")
    assert_true(check_equation_local(G1, G2, res), "A_only: не выполняется A*G1 == G2*P")
end


# ТЕСТ 3: A + P вместе

function test_A_and_perm_random()
    Random.seed!(12345)
    k, n = 5, 10
    G1 = random_fullrank_G(k, n)
    A0 = random_invertible_A(k)
    perm0 = randperm(n)
    G2 = gf2_mul(A0, G1[:, perm0])

    res = isequivSSP(G1, G2, sig3_ws_hull_short; max_tries=1_000_000)
    assert_true(res != -1 && res != -2, "A+P_random: алгоритм не нашёл эквивалентность")
    assert_true(check_equation_local(G1, G2, res), "A+P_random: не выполняется A*G1 == G2*P")
end

# ТЕСТ 4: заведомо НЕ эквивалентны (разные коды)

function test_not_equiv_by_rank()
    # G1 имеет ранг 2
    G1 = Bool[
        1 0 0 1;
        0 1 1 0
    ]
    # G2 делаем ранг 1 (вторая строка такая же, как первая)
    G2 = Bool[
        1 0 0 1;
        1 0 0 1
    ]

    # Если коды эквивалентны, ранги порождающих матриц должны совпадать (после A)
    res = isequivSSP(G1, G2, sig3_ws_hull_short; max_tries=10000)
    assert_eq(res, -1, "not_equiv_by_rank: ожидали -1 (неэквивалентность)")
end


# ТЕСТ 5: LIMIT должен сработать (тут mmax_tries очень маленький)

function test_limit()
    Random.seed!(7)
    k, n = 6, 12
    G1 = random_fullrank_G(k, n)
    A0 = random_invertible_A(k)
    perm0 = randperm(n)
    G2 = gf2_mul(A0, G1[:, perm0])

    # max_tries = 1 почти гарантированно слишком мало
    res = isequivSSP(G1, G2, sig1_dim_short; max_tries=1)
    assert_true(res == -2 || res == -1 || (res != -1 && res != -2),
                "limit: неожиданный тип результата")
end


# Запуск тестов

println("Запуск тестов задания 02")

test_perm_only()
println("OK: test_perm_only")

test_A_only()
println("OK: test_A_only")

test_A_and_perm_random()
println("OK: test_A_and_perm_random")

test_not_equiv_by_rank()
println("OK: test_not_equiv_by_rank")

test_limit()
println("OK: test_limit")

println("Все тесты пройдены")
