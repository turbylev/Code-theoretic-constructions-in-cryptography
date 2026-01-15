using Random
include("isequivSSP.jl")

println("Задание 02 (SSA)")

# Введите размерность порождающей матрицы (n x k)
n = 25
k = 12
#k = n ÷ 2

# Введите число экспериментов для каждой сигнатуры для вычисления среднего времени работы
REPS = 1

# Введите лимит перебора, то есть максимальное число полных perm, для которых пробуем найти A


println("Параметры: k=$k, n=$n, прогонов=$REPS")

# Проверка равенства A*G1 == G2*P
function check_equation(G1, G2, res)
    if res == -1 || res == -2
        return false
    end
    A, perm = res
    P = perm_to_P(perm)
    left  = gf2_mul(toBoolMat(A), toBoolMat(G1))
    right = gf2_mul(toBoolMat(G2), P)
    return left == right
end

# Функция генерации рандомной матрицы и эквивалентной пары
function random_G(k, n)
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

function generate_pair(k, n)
    G1 = random_G(k, n)  # <-- исправлено: было random_fullrank_G
    A  = random_invertible_A(k)
    perm = randperm(n)
    G2 = gf2_mul(A, toBoolMat(G1)[:, perm])
    return G1, G2
end

G1, G2 = generate_pair(k, n)

println("\nСгенерировали G1 и эквивалентную G2.")
println("G1:"); print01(G1)
println("G2:"); print01(G2)


# функция серии прогонов для одной сигнатуры
function run_series(sigfun, name; reps=20, max_tries=50000)
    println("\n Эквивалентность ($name)")
    println("Лимит перебора max_tries = $max_tries на каждый прогон")

    total = 0.0
    done = 0
    last_res = -1
    stopped_by_limit = false

    for t in 1:reps
        elapsed = @elapsed begin
            last_res = isequivSSP(G1, G2, sigfun; max_tries=max_tries)
        end

        if last_res == -2
            println("Ошибка: прогон $t остановлен по max_tries.")
            stopped_by_limit = true
            break
        end

        total += elapsed
        done += 1
    end

    if done == 0
        print_result(-2)
        return Inf, -2
    end

    avg = total / done

    print_result(last_res)

    if last_res != -1 && last_res != -2
        ok = check_equation(G1, G2, last_res)
        println("Проверка A*G1 == G2*P: ", ok ? "успешно" : "ошибка")
    end

    println("Выполнено прогонов: $done / $reps")
    println("Среднее время: ", round(avg, digits=6), " сек")
    if stopped_by_limit
        println("Warning: серия остановлена из-за лимита перебора.")
    end

    return avg, last_res
end

t1, _ = run_series(sig1_dim_short,      "sig1 (dim укороченного кода)"; reps=REPS, max_tries=50000)
t2, _ = run_series(sig2_dim_hull_short, "sig2 (dim hull укороченного)"; reps=REPS, max_tries=200000)
t3, _ = run_series(sig3_ws_hull_short,  "sig3 (спектр весов hull укороченного)"; reps=REPS, max_tries=50000000)

println("\nИтог по среднему времени для матриц размера n = ", n, " k = ", k)
println("sig1: ", isfinite(t1) ? string(round(t1, digits=6), " сек") : "LIMIT")
println("sig2: ", isfinite(t2) ? string(round(t2, digits=6), " сек") : "LIMIT")
println("sig3: ", isfinite(t3) ? string(round(t3, digits=6), " сек") : "LIMIT")
