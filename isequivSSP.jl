# isequivSSP.jl
# Задание 02: “Ой, Кеша, они же эквивалентны!”
#
# Проверяем эквивалентность порождающих матриц G1 и G2 (k×n):
#   существует A ∈ GL(k,2) и перестановочная матрица P (n×n), такие что
#       A * G1 == G2 * P
#
# Здесь P представляем через вектор perm (перестановка столбцов):
#       G2[:, perm] == G2 * P
#
# Главная функция:
#   isequivSSP(G1, G2, sigfun; max_tries=...) -> -1 / -2 / (A, perm)
#
# Возвращает:
#   -1  : не найдено (коды не эквивалентны ИЛИ сигнатура не помогла)
#   -2  : достигнут лимит перебора перестановок (max_tries)
#   (A, perm): найдено преобразование

# Всякие доп функции

toBoolMat(G) = Bool.(G .% 2) # приведение к GF(2) в виде Bool
wt(v) = count(identity, v)  # вес (число единиц)

function print01(M)
    A = Int.(M)
    if ndims(A) == 1
        println(join(A, " "))
        return
    end
    for i in 1:size(A, 1)
        println(join(A[i, :], " "))
    end
end

# строит перестановочную матрицу P (n×n) так, что G[:,perm] == G * P
function perm_to_P(perm)
    n = length(perm)
    P = falses(n, n)
    for j in 1:n
        P[perm[j], j] = true
    end
    return P
end

function print_result(res)
    if res == -1
        println("Результат: коды не эквивалентны или сигнатура не дала ответа")
        return
    end
    if res == -2
        println("Результат: достигнут лимит перебора перестановок (max_tries).")
        return
    end
    A, perm = res
    println("Результат: коды эквивалентны.")
    println("Матрица A:")
    print01(A)
    println("Перестановка perm = ", perm)
    println("Перестановочная матрица P:")
    P = perm_to_P(perm)
    print01(P)
end

# функция умножение матриц над GF(2)

function gf2_mul(A, B)
    r, m = size(A)
    m2, c = size(B)
    m == m2 || error("gf2_mul: размерности матриц не совпадают")
    C = falses(r, c)
    @inbounds for i in 1:r
        for k in 1:m
            if A[i, k]
                for j in 1:c
                    C[i, j] ⊻= B[k, j]
                end
            end
        end
    end
    return C
end

# Гаусс / RREF / ранг над GF(2)

function rref_gf2(A)
    R = copy(A)
    rows, cols = size(R)
    pivcols = Int[]
    r = 1
    @inbounds for c in 1:cols
        r > rows && break

        piv = 0
        for rr in r:rows
            if R[rr, c]
                piv = rr
                break
            end
        end
        piv == 0 && continue

        if piv != r
            R[r, :], R[piv, :] = R[piv, :], R[r, :]
        end
        push!(pivcols, c)

        for rr in 1:rows
            if rr != r && R[rr, c]
                for j in c:cols
                    R[rr, j] ⊻= R[r, j]
                end
            end
        end

        r += 1
    end
    return R, pivcols, length(pivcols)
end

rank_gf2(A) = last(rref_gf2(A))

# функция решения M*x=b над GF(2)
function solve_gf2(M, b)
    r, n = size(M)
    length(b) == r || error("solve_gf2: несовместимые размеры M и b")

    Aug = falses(r, n + 1)
    Aug[:, 1:n] .= M
    Aug[:, n + 1] .= b

    R, pivcols, _ = rref_gf2(Aug)

    @inbounds for i in 1:r
        allzero = true
        for j in 1:n
            if R[i, j]
                allzero = false
                break
            end
        end
        if allzero && R[i, n + 1]
            return nothing
        end
    end

    x = falses(n)
    prow = 1
    @inbounds for pcol in pivcols
        pcol > n && continue
        x[pcol] = R[prow, n + 1]
        prow += 1
        prow > r && break
    end
    return x
end

# базис ядра: решения A*x=0 (x — столбец). Возвращает матрицу B (d×n), где здесь строки - это базисные векторы.
function nullspace_basis(A)
    R, pivcols, _ = rref_gf2(A)
    _, n = size(A)

    pivset = Set(pivcols)
    freecols = [j for j in 1:n if !(j in pivset)]
    d = length(freecols)

    B = falses(d, n)
    d == 0 && return B

    @inbounds for (idx, fcol) in enumerate(freecols)
        x = falses(n)
        x[fcol] = true
        for (prow, pcol) in enumerate(pivcols)
            s = false
            for fc in freecols
                if R[prow, fc] && x[fc]
                    s ⊻= true
                end
            end
            x[pcol] = s
        end
        B[idx, :] .= x
    end
    return B
end

is_invertible_gf2(A) = (size(A,1) == size(A,2) && rank_gf2(A) == size(A,1))


# функция укорочение по позиции i

function shortened_generator(G, i)
    k, n = size(G)
    (1 <= i <= n) || error("shortened_generator: i не вхдит в диапазон координат]")
    g = G[:, i]

    if !any(g)
        B = falses(k, k)
        for t in 1:k
            B[t, t] = true
        end
        return gf2_mul(B, G)
    end

    p = findfirst(==(true), g)
    free = [j for j in 1:k if j != p]

    B = falses(k - 1, k)
    @inbounds for (row, j) in enumerate(free)
        B[row, j] = true
        if g[j]
            B[row, p] = true
        end
    end
    return gf2_mul(B, G)
end

# Hull
function hull_generator(G)
    S = gf2_mul(G, transpose(G))
    B = nullspace_basis(S)
    return gf2_mul(B, G)
end

# спектр весов
function weight_spectrum_from_generator(H)
    h, n = size(H)
    D = zeros(Int, n + 1)
    if h == 0
        D[1] = 1
        return D
    end
    maxT = (UInt(1) << h) - UInt(1)
    for t in UInt(0):maxT
        v = falses(n)
        @inbounds for r in 1:h
            if ((t >> (r - 1)) & UInt(1)) == UInt(1)
                for j in 1:n
                    v[j] ⊻= H[r, j]
                end
            end
        end
        D[wt(v) + 1] += 1
    end
    return D
end


# функции для сигнатур

sig1_dim_short(G, i) = rank_gf2(shortened_generator(G, i))

function sig2_dim_hull_short(G, i)
    Gs = shortened_generator(G, i)
    Hs = hull_generator(Gs)
    return rank_gf2(Hs)
end

function sig3_ws_hull_short(G, i)
    Gs = shortened_generator(G, i)
    Hs = hull_generator(Gs)
    return Tuple(weight_spectrum_from_generator(Hs))
end


# разбиение позиций по сигнатуре

function signature_classes(G, sigfun)
    _, n = size(G)
    cls = Dict{Any, Vector{Int}}()
    for i in 1:n
        s = sigfun(G, i)
        if haskey(cls, s)
            push!(cls[s], i)
        else
            cls[s] = [i]
        end
    end
    return cls
end


# функция нахождения A при фиксированной перестановке столбцов G2
function solve_for_A(G1, G2perm)
    k, n = size(G1)
    size(G2perm) == (k, n) || return nothing

    M = transpose(G1)        # n×k
    A = falses(k, k)

    @inbounds for r in 1:k
        b = vec(transpose(G2perm[r, :]))
        x = solve_gf2(M, b)
        x === nothing && return nothing
        A[r, :] .= x
    end

    is_invertible_gf2(A) || return nothing
    gf2_mul(A, G1) == G2perm || return nothing
    return A
end


# перебор перестановок
function permute_first!(cb, arr, l)
    if l > length(arr)
        return cb(arr)
    end
    for i in l:length(arr)
        arr[l], arr[i] = arr[i], arr[l]
        ans = permute_first!(cb, arr, l + 1)
        arr[l], arr[i] = arr[i], arr[l]
        ans !== nothing && return ans
    end
    return nothing
end

# перебор перестановок (НО с возможностью мгновенно остановиться по флагу)
function permute_first_stop!(cb, arr, l, should_stop)
    should_stop() && return nothing

    if l > length(arr)
        return cb(arr)
    end

    for i in l:length(arr)
        should_stop() && return nothing

        arr[l], arr[i] = arr[i], arr[l]
        ans = permute_first_stop!(cb, arr, l + 1, should_stop)
        arr[l], arr[i] = arr[i], arr[l]

        ans !== nothing && return ans
    end
    return nothing
end


# основная функция из задания isequivSSP

function isequivSSP(G1_in, G2_in, sigfun; max_tries::Int=typemax(Int))
    G1 = toBoolMat(G1_in)
    G2 = toBoolMat(G2_in)

    k1, n1 = size(G1)
    k2, n2 = size(G2)
    (k1 == k2 && n1 == n2) || return -1
    n = n1

    cls1 = signature_classes(G1, sigfun)
    cls2 = signature_classes(G2, sigfun)

    if length(cls1) != length(cls2)
        return -1
    end
    for (s, idx1) in cls1
        if !haskey(cls2, s) || length(cls2[s]) != length(idx1)
            return -1
        end
    end

    sigs = collect(keys(cls1))
    sort!(sigs, by = s -> -length(cls1[s]))

    perm = fill(0, n)

    tries = 0
    limit_reached = false

    function backtrack(t)
        if limit_reached
            return nothing
        end

        if t > length(sigs)
            tries += 1
            if tries > max_tries
                limit_reached = true
                return nothing
            end

            G2perm = G2[:, perm]
            A = solve_for_A(G1, G2perm)
            return A === nothing ? nothing : (A, copy(perm))
        end

        s = sigs[t]
        idx1 = cls1[s]
        idx2 = cls2[s]
        arr = copy(idx2)

        should_stop() = limit_reached

        return permute_first_stop!((current) -> begin
            if limit_reached
                return nothing
            end
            @inbounds for j in 1:length(idx1)
                perm[idx1[j]] = current[j]
            end
            return backtrack(t + 1)
        end, arr, 1, should_stop)
    end

    ans = backtrack(1)
    if ans === nothing
        return limit_reached ? -2 : -1
    end
    return ans
end
