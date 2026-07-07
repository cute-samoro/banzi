#set page(
  paper: "a4",
  margin: (top: 1.2cm, bottom: 1.4cm, x: 1.25cm),
  numbering: "1",
  number-align: center + bottom,
)

#set text(size: 9pt)
#set par(leading: 0.52em)
#set heading(numbering: "1.")
#set raw(theme: auto, tab-size: 2)

#show raw: set text(
  font: "Cascadia Code",
  size: 0.92em,
  fill: rgb("#111827"),
)

#show raw.where(block: false): box.with(
  fill: rgb("#f3f4f6"),
  inset: (x: 0.35em, y: 0.05em),
  radius: 3pt,
)

#show raw.where(block: true): block.with(
  fill: rgb("#f8fafc"),
  stroke: 0.45pt + rgb("#e5e7eb"),
  inset: (x: 10pt, y: 8pt),
  radius: 6pt,
)

#align(center)[
  #text(16pt, weight: "bold")[算法板子]
]

#outline(depth: 2)
#pagebreak()
#set page(columns: 1)

= #text("NTT")
```cpp
const int MOD = 998244353, G = 3; // 常见 NTT 模数及其原根
//power是ksm
void ntt(vector<int>& a, bool inv) {
    int n = a.size();
    for (int i = 1, j = 0; i<n; i++) {
        int bit = n>>1;
        for (; j&bit; bit>>= 1) j ^= bit;
        j ^= bit;
        if (i<j) swap(a[i], a[j]);
    }
    for (int len = 2; len<= n; len<<= 1) {
        ll wlen = power(G, (MOD - 1) / len);
        if (inv) wlen = power(wlen, MOD - 2);
        for (int i = 0; i<n; i += len) {
            ll w = 1;
            for (int j = 0; j<len / 2; j++) {
                int u = a[i + j], v = (int)(1LL * a[i + j + len / 2] * w % MOD);
                a[i + j] = (u + v) % MOD;
                a[i + j + len / 2] = (u - v + MOD) % MOD;
                w = w * wlen % MOD;
            }
        }
    }
    if (inv) {
        ll n_inv = power(n, MOD - 2);
        for (int&x : a) x = (int)(1LL * x * n_inv % MOD);
    }
}
vector<int>multiply(vector<int>a, vector<int>b) {
    int n = 1, total = a.size() + b.size() - 1;
    while (n<total) n<<= 1;
    a.resize(n); b.resize(n);
    ntt(a, 0); ntt(b, 0);
    for (int i = 0; i<n; i++) a[i] = 1LL * a[i] * b[i] % MOD;
    ntt(a, 1);
    a.resize(total);
    //求循环卷积return res
    /*vector<int>res(L, 0);
    for (int i = 0; i<n; i++) {
        res[i % L] = (res[i % L] + a[i]) % MOD;
    }*/
    return a;}
```

= #text("凸包（andrew算法）")

```cpp
struct Point {
    long long x, y;
    // 排序：先按 x 升序，再按 y 升序
    bool operator<(const Point& t) const {
        return x == t.x ? y<t.y : x<t.x;
    }
};
// 叉积：(B-A) x (C-A)。>0 表示 C 在 AB 左侧（逆时针旋转）
long long cross(Point a, Point b, Point c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}
vector<Point>getConvexHull(vector<Point>&p) {
    int n = p.size(), k = 0;
    if (n <= 2) return p;
    vector<Point>res(2 * n);
    sort(p.begin(), p.end());
    // 求下凸壳
    for (int i = 0; i<n; i++) {
        while (k>= 2&&cross(res[k - 2], res[k - 1], p[i])<= 0) k--;
        res[k++] = p[i];
    }
    // 求上凸壳
    for (int i = n - 2, t = k + 1; i>= 0; i--) {
        while (k>= t&&cross(res[k - 2], res[k - 1], p[i])<= 0) k--;
        res[k++] = p[i];
    }
    res.resize(k - 1); // 最后一个点和第一个点重合，删掉
    return res;
}
```

= #text("Exgcd")

```cpp
i64 exgcd(i64 a, i64 b, i64&x, i64&y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i64 g = exgcd(b, a % b, y, x);
    y -= a / b * x;
    return g;
}
```

= #text("欧拉路径")

```cpp
const int MAXN = 100005;
struct Edge {
    int to, id;
};
vector<Edge>adj[MAXN];
int head[MAXN]; // 当前弧优化：记录每个点处理到了第几条边
bool used[MAXN * 2]; // 标记边是否被访问过（处理无向图或多重边）
vector<int> path; // 存储最终路径（点序列或边序列）

void dfs(int u) {
    // 遍历从 u 出发的边，head[u] 保证每条边只被扫描一次
    for (int &i = head[u]; i<adj[u].size(); ) {
        Edge e = adj[u][i++]; // 注意这里必须先 i++ 再递归
        if (used[e.id]) continue;
        used[e.id] = true; // 标记该边已用
        dfs(e.to);
    }
    path.push_back(u); // 回溯时加入序列
}

// 得到路径：reverse(path.begin(), path.end());
```

= #text("马拉车")

```cpp
// d[i] 表示以 i 为中心的最长回文半径
// 转换后的字符串长度为 2n+1，d[i]-1 即为原字符串中对应回文串的长度
vector<int>manacher(string s) {
    string t = "#";
    for (char c : s) t += c, t += "#";
    int n = t.size();
    vector<int>d(n);
    // r 为当前已知回文串覆盖的最右边界，l 为该边界对应的中心
    for (int i = 0, l = 0, r = -1; i<n; i++) {
        int k = (i>r) ? 1 : min(d[l + r - i], r - i + 1);
        while (0<= i - k&&i + k<n&&t[i - k] == t[i + k]) k++;
        d[i] = k--;
        if (i + k>r) l = i - k, r = i + k;
    }
    return d;
}
```

= #text("欧拉筛")

```cpp
const int MAXN = 1000005;
int primes[MAXN], cnt;
bool is_not_prime[MAXN];

void sieve(int n) {
    for (int i = 2; i<= n; i++) {
        if (!is_not_prime[i]) primes[cnt++] = i;
        for (int j = 0; j<cnt&&i * primes[j]<= n; j++) {
            is_not_prime[i * primes[j]] = true;
            if (i % primes[j] == 0) break; // 核心：保证每个合数只被最小质因子筛掉
        }
    }
}
```

= #text("拉格朗日插值（通用版，o(n²)）")

```cpp
// n 个点 (x[i], y[i])，求 f(k)
ll lagrange(int n, ll k, vector<ll>&x, vector<ll>&y) {
    ll ans = 0;
    for (int i = 0; i<n; i++) {
        ll up = y[i], down = 1;
        for (int j = 0; j<n; j++) {
            if (i == j) continue;
            up = up * (k - x[j] % MOD + MOD) % MOD;
            down = down * (x[i] - x[j] % MOD + MOD) % MOD;
        }
        ans = (ans + up * inv(down)) % MOD;
    }
    return ans;
}
```

= #text("拉格朗日插值（连续整数）")

```cpp
ll pre[MAXN], suf[MAXN], fac[MAXN], inv_fac[MAXN];
// x_i 为 1, 2, ..., n 的情况，求 f(k)
ll lagrange_linear(int n, ll k, vector<ll>&y) {
    if (k >= 1 && k <= n) return y[k - 1];

    pre[0] = suf[n + 1] = 1;
    for (int i = 1; i<= n; i++) pre[i] = pre[i - 1] * (k - i % MOD + MOD) % MOD;
    for (int i = n; i>= 1; i--) suf[i] = suf[i + 1] * (k - i % MOD + MOD) % MOD;

    fac[0] = 1;
    for (int i = 1; i<= n; i++) fac[i] = fac[i - 1] * i % MOD;
    inv_fac[n] = inv(fac[n]);
    for (int i = n - 1; i>= 0; i--) inv_fac[i] = inv_fac[i + 1] * (i + 1) % MOD;

    ll ans = 0;
    for (int i = 1; i<= n; i++) {
        ll up = pre[i - 1] * suf[i + 1] % MOD * y[i - 1] % MOD;
        ll down = inv_fac[i - 1] * inv_fac[n - i] % MOD;
        // 分母符号取决于 (n-i) 的奇偶性
        if ((n - i) % 2) ans = (ans - up * down % MOD + MOD) % MOD;
        else ans = (ans + up * down % MOD) % MOD;
    }
    return ans;
}
```

= #text("Bsgs（gcd(a,p)==1）")

```cpp
ll bsgs(ll a, ll b, ll p) {
    if (1 % p == b % p) return 0;
    unordered_map<ll, ll>hash;
    ll m = ceil(sqrt(p));

    // Baby-step: 预处理 a^j * b 并存入哈希表
    ll cur = b % p;
    for (int j = 0; j<m; j++) {
        hash[cur] = j;
        cur = cur * a % p;
    }

    // Giant-step: 计算 am = a^m，然后枚举 i 检查 am^i
    ll am = 1;
    for (int i = 0; i<m; i++) am = am * a % p;

    cur = am;
    for (int i = 1; i<= m; i++) {
        if (hash.count(cur)) return i * m - hash[cur];
        cur = cur * am % p;
    }
    return -1; // 无解
}
```

= #text("Excrt（扩展中国剩余定理）")

```cpp
// 这里的乘法取模是为了防止 (a * b) % m 溢出，如果 a, b 较小可直接使用 a * b % m
ll qmul(ll a, ll b, ll m) {
    ll res = 0;
    while (b) {
        if (b&1) res = (res + a) % m;
        a = (a + a) % m;
        b>>= 1;
    }
    return res;
}
// n 个方程：x = r[i] (mod m[i])
ll excrt(vector<ll>&m, vector<ll>&r) {
    ll M = m[0], R = r[0]; // R 为当前的解，M 为当前的最小公倍数
    for (int i = 1; i<m.size(); i++) {
        ll x, y;
        ll d = exgcd(M, m[i], x, y);
        if ((r[i] - R) % d != 0) return -1; // 无解

        // 求解 M*x = r[i] - R (mod m[i])
        ll mod = m[i] / d;
        x = qmul((r[i] - R) / d, x, mod); // 这里的乘法注意溢出
        if (x<0) x += mod;

        R += x * M;
        M = M / d * m[i];
        R = (R % M + M) % M;
    }
    return R;
}
```

= #text("Odt")

```cpp
struct Node {
    int l, r;
    mutable ll v; // mutable 允许在 set 中直接修改 v
    Node(int L, int R = -1, ll V = 0) : l(L), r(R), v(V) {}
    bool operator<(const Node& t) const { return l<t.l; }
};

set<Node>odt;

// 分裂操作：将包含 pos 的区间 [l, r] 分裂为 [l, pos-1] 和 [pos, r]
// 返回后一个区间的迭代器
auto split(int pos) {
    auto it = odt.lower_bound(Node(pos));
    if (it != odt.end() && it->l == pos) return it;
    --it;
    int l = it->l, r = it->r;
    ll v = it->v;
    odt.erase(it);
    odt.insert(Node(l, pos - 1, v));
    return odt.insert(Node(pos, r, v)).first;
}

// 核心：区间赋值（推平操作）
void assign(int l, int r, ll v) {
    auto itR = split(r + 1), itL = split(l); // 顺序必须先 R 后 L
    odt.erase(itL, itR);
    odt.insert(Node(l, r, v));
}

// 其他区间操作模板（以区间加为例）
void add(int l, int r, ll v) {
    auto itR = split(r + 1), itL = split(l);
    for (; itL != itR; ++itL) itL->v += v;
}
```

= #text("前缀数组（来自oiwiki）")

```cpp
vector<int>prefix_function(string s) {
  int n = (int)s.length();
  vector<int>pi(n);
  for (int i = 1; i<n; i++) {
    int j = pi[i - 1];
    while (j>0&&s[i] != s[j]) j = pi[j - 1];
    if (s[i] == s[j]) j++;
    pi[i] = j;
  }
  return pi;
}
```

= #text("kmp（来自oiwiki）")

```cpp
vector<int>find_occurrences(string text, string pattern) {
  string cur = pattern + '#' + text;
  int sz1 = text.size(), sz2 = pattern.size();
  vector<int>v;
  vector<int>lps = prefix_function(cur);
  for (int i = sz2 + 1; i<= sz1 + sz2; i++) {
    if (lps[i] == sz2) v.push_back(i - 2 * sz2);
  }
  return v;
}
```

= #text("mt19937")

```cpp
#include<bits/stdc++.h>
int main() {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<>dis(1, 100);
    //std::uniform_real_distribution<double>dis(0.0, 1.0);
    for (int i = 0; i<5; ++i) {
        std::cout << dis(gen) << " ";
    }

    return 0;
}
```

= #text("二分图最大匹配（匈牙利算法，稠密图且点数大于1e4不适用）")

```cpp
const int MAXN = 505; // 左侧集合点的最大数量
vector<int> adj[MAXN]; // 邻接表，只存从左往右的边
int match[MAXN];      // match[y] = x 表示右侧点 y 匹配了左侧点 x
bool vis[MAXN];       // 标记右侧点在单次 DFS 中是否被访问过
bool dfs(int u) {
    for (int v : adj[u]) {
        if (!vis[v]) {
            vis[v] = true;
            // 如果右侧点没有匹配，或者原匹配点可以找到新的增广路
            if (match[v] == -1 || dfs(match[v])) {
                match[v] = u;
                return true;
            }
        }
    }
    return false;
}
int solve(int n_left) {
    int ans = 0;
    memset(match, -1, sizeof(match));
    for (int i = 1; i<= n_left; i++) {
        memset(vis, false, sizeof(vis)); // 每次都要重置
        if (dfs(i)) ans++;
    }
    return ans;
}
```

= #text("线性基")

```cpp
typedef long long ll;
const int MAXL = 60;
struct LinearBasis {
    ll a[MAXL + 1];

    LinearBasis() {
        memset(a, 0, sizeof(a));
    }

    // 简洁插入（不消元，只保证线性无关）
    void insert(ll x) {
        for (int i = MAXL; i>= 0; --i) {
            if (!(x >> i & 1)) continue;
            if (!a[i]) { a[i] = x; return; }
            x ^= a[i];
        }
    }

    // 查询异或最大值（贪心即可，不需要消元）
    ll queryMax() {
        ll res = 0;
        for (int i = MAXL; i>= 0; --i)
            if ((res ^ a[i])>res) res ^= a[i];
        return res;
    }

    // 合并另一个线性基
    void merge(const LinearBasis&other) {
        for (int i = 0; i<= MAXL; ++i)
            if (other.a[i]) insert(other.a[i]);
    }
};
int main() {
    int n;
    cin>>n;
    LinearBasis lb;
    for (int i = 0; i<n; ++i) {
        ll x;
        cin>>x;
        lb.insert(x);
    }
    cout << lb.queryMax() << "\n";
    return 0;}
```

= #text("强连通分量")

```cpp
const int N = 1e5 + 5;
vector<int>g[N];
int dfn[N], low[N], timer;
bool inStk[N];
vector<int> stk;                     // 模拟栈
vector<vector<int>> scc;            // 收集所有强连通分量

void tarjan(int u) {
    dfn[u] = low[u] = ++timer;
    stk.push_back(u);
    inStk[u] = true;
    for (int v : g[u]) {
        if (!dfn[v]) {
            tarjan(v);
            low[u] = min(low[u], low[v]);
        } else if (inStk[v]) {
            low[u] = min(low[u], dfn[v]);
        }
    }
    if (dfn[u] == low[u]) {          // u 是 SCC 的根
        vector<int>comp;
        while (true) {
            int v = stk.back();
            stk.pop_back();
            inStk[v] = false;
            comp.push_back(v);
            if (v == u) break;
        }
        scc.push_back(comp); }}        // 收集当前分量
int main() {
    read(n,m,g);//g为edge
    for (int i = 1; i<= n; ++i)
        if (!dfn[i]) tarjan(i);
    // 输出所有强连通分量
    for (int i = 0; i<(int)scc.size(); ++i) {
        cout << "SCC " << i + 1 << ": ";
        for (int x : scc[i]) cout << x << " ";
        cout << "\n";
    }
    return 0;}
```

= #text("polya定理")
$|C/G| = frac(1,|G|) sum_(g in G) m^c(g)$
- #text("置换群")
- #text("群的阶（置换个数）")
- #text("所有着色方案（共 m^n 种）")
- #text("轨道集合（本质不同着色）")
- #text("可用颜色数")
- #text("置换 g 的循环个数")
- #text("在 g 下保持不变的着色数")

$|X/G| = frac(1, |G|) sum_(g in G) |X^g|$

#text("补充：burnside引理")
- #text("被作用的集合")
- #text("轨道集")
- #text("群阶")
- #text("g 作用下的不动点集")

= #text("求割点")

```cpp
const int N = 1e5 + 5;
vector<int>G[N];
int dfn[N], low[N], timer;
bool cut[N];           // cut[u] = true 表示 u 是割点

void tarjan(int u, int fa) {
    dfn[u] = low[u] = ++timer;
    int child = 0;                     // 子树个数（仅对根有用）
    for (int v : G[u]) {
        if (v == fa) continue;
        if (!dfn[v]) {
            child++;
            tarjan(v, u);
            low[u] = min(low[u], low[v]);
            if (low[v] >= dfn[u]) cut[u] = true;  // 割点判定
        } else {
            low[u] = min(low[u], dfn[v]);
        }
    }
    if (fa == -1) cut[u] = (child >= 2);  // 根节点特判
}
// 调用：
// for (int i = 1; i<= n; i++)
//     if (!dfn[i]) tarjan(i, -1);
```

= #text("自定义hash")

```cpp
struct Person {
    std::string name;
    int age;
};
// 必须在 std 中打开 namespace
namespace std {
    template<>
    struct hash<Person>{
        size_t operator()(const Person&p) const noexcept {
            // 组合多个成员的哈希值（常用方法见后文）
            size_t h1 = hash<string>{}(p.name);
            size_t h2 = hash<int>{}(p.age);
            return h1 ^ (h2 << 1);   // 简单组合，更稳健的方法见下
        }
    };
}
```

```cpp
namespace std {
    template<typename T1, typename T2>
    struct hash<pair<T1, T2>>{
        size_t operator()(const pair<T1, T2>&p) const {
            auto h1 = hash<T1>{}(p.first);
            auto h2 = hash<T2>{}(p.second);
            // 组合方法...
        }
    };
}
```

```cpp
struct Point {
    int x, y;
};
namespace std {
    template<>
    struct hash<Point>{
        size_t operator()(const Point&p) const {
            size_t seed = 0;
            hash<int>hasher;
            // 对每个字段依次调用 hash_combine
            hash_combine(seed, hasher(p.x));
            hash_combine(seed, hasher(p.y));
            return seed;
        }
    };
}
inline void hash_combine(size_t& seed, size_t hash) {
    seed ^= hash + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}
```

= #text("hall定理")
- #text("设二分图 G = (L, R, E)，L 为左部，R 为右部。")
- #text("1) 存在匹配覆盖 L 中所有顶点 当且仅当 对任意子集 S ⊆ L，其邻域 N(S) 满足：")
- #text("|N(S)| ≥ |S|")
- #text("2) 最大匹配大小 = |L| - max_{( |S| - |N(S)| ) | S ⊆ L}")
- #text("（若该值 < |L|，则无法完全覆盖 L）")

= #text("Zobrist Hashing（xor hashing）")

#text("赋值：每个值 v → 随机 64 位整数 hash[v]")

#text("前缀：pre[i] = pre[i-1] ^ hash[A[i]]")

#text("比较：子数组 XOR = pre[r] ^ pre[l-1]，与预期哈希比较")

#text("多测放全局开mt19937")

```cpp
//mt19937 rng32(...);    // 生成 uint32_t
//mt19937_64 rng64(...); // 生成 uint64_t
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
hash[i] = rng();
```

= #text("Lemma (Bertrand's postulate)") 
```text
For each positive integer xx, there is a prime pp inside the interval [x,2x].
```

= #text("求多个数的欧拉函数")

```cpp
vector<int> cal_euler(int x) {
    vector<int> re(x + 1);
    vector<int> prime;
    vector<int> not_prime(x + 1, 0);
    for (int i = 2; i <= x; i++) {
        if (!not_prime[i]) {
            prime.emplace_back(i);
            re[i] = i - 1;
        }
        for (auto j : prime) {
            if (i * j > x) break;
            not_prime[i * j] = 1;
            if (i % j == 0) {
                re[i * j] = re[i] * j;
            }
            else {
                re[i * j] = re[i] * re[j];
            }
            if (i % j == 0) break;
        }
    }
    return re;
}
```

= #text("最大流（dinic算法）")

```cpp
struct Dinic {
    struct Edge { int to, rev; long long cap; };
    vector<vector<Edge>> g;
    vector<int> lev, iter;
    void init(int n) { g.assign(n, {}); }
    void add_edge(int u, int v, long long cap) {
        g[u].push_back({v, (int)g[v].size(), cap});
        g[v].push_back({u, (int)g[u].size() - 1, 0});
    }
    void bfs(int s) {
        fill(lev.begin(), lev.end(), -1);
        queue<int> q;
        lev[s] = 0; q.push(s);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (auto &e : g[u]) {
                if (e.cap > 0 && lev[e.to] == -1) {
                    lev[e.to] = lev[u] + 1;
                    q.push(e.to);
                }
            }
        }
    }
long long dfs(int u, int t, long long f) {
      if (u == t) return f;
      for (int &i = iter[u]; i < (int)g[u].size(); ++i) {
          auto &e = g[u][i];
          if (e.cap > 0 && lev[e.to] == lev[u] + 1) {
              long long d = dfs(e.to, t, min(f, e.cap));
              if (d > 0) {
                  e.cap -= d;
                  g[e.to][e.rev].cap += d;
                  return d;
              }
          }
      }
      return 0;
  }
long long max_flow(int s, int t) {
        long long flow = 0;
        lev.resize(g.size()); iter.resize(g.size());
        while (true) {
            bfs(s);
            if (lev[t] == -1) break;
            fill(iter.begin(), iter.end(), 0);
            long long f;
            while ((f = dfs(s, t, 1e18)) > 0) flow += f;
        }
        return flow;
    }
};
最小割：最大流 = 最小割。求完最大流后，从源点沿剩余容量 >0 的边 DFS，能到达的点集就是 S，其余为 T。
```

= #text("有源汇最大流")

```cpp
struct BoundFlow {
    Dinic dinic;
    vector<long long> in;
    int S, T; // 超级源汇
    void init(int n) {
        dinic.init(n + 2);
        in.assign(n + 2, 0);
        S = n; T = n + 1;
    }
    // 添加有上下界的边 (u->v, l, r)
    void add_edge(int u, int v, long long l, long long r) {
        in[u] -= l; in[v] += l;
        dinic.add_edge(u, v, r - l);
    }
    bool feasible() {
        long long total = 0;
        for (int i = 0; i < S; ++i) { // S 是原图中最后一个点下标+1 ？
            if (in[i] > 0) {
                dinic.add_edge(S, i, in[i]);
                total += in[i];
            } else if (in[i] < 0) {
                dinic.add_edge(i, T, -in[i]);
            }
        }
        return dinic.max_flow(S, T) == total;
    }
// 有源汇 s,t 的最大流（先调用 add_edge 建图再调这个）
    long long max_flow_with_bound(int s, int t) {
        add_edge(t, s, 0, 1e18); // 加下限0的边
        if (!feasible()) return -1; // 无解
        long long flow = dinic.g[s].back().cap; // t->s 反向边容量 = 可行流中 s->t 的流量
        // 删掉 t->s 附加边
        dinic.g[t].pop_back(); dinic.g[s].pop_back();
        return flow + dinic.max_flow(s, t);
    }
};
上面 feasible() 中的 S 是 S 的点编号，初始化时 S = n，但循环里 i < S 跑到了 n-1，即所有原图点（0~n-1）。初始化时传入点数即可。
无源汇可行流
每个点计算 in[v] = 所有入边下限和 - 出边下限和。
建超级源 S 和超级汇 T。
若 in[v] > 0，连 S -> v 容量 in[v]；若 in[v] < 0，连 v -> T 容量 -in[v]。
原边连容量 r - l。
跑最大流，若从 S 出发的所有边满流，则有解
```

= #text("最小费用最大流（spfa）")

```cpp
struct MCMF {
    struct Edge { int to, rev; long long cap, cost; };
    vector<vector<Edge>> g;
    vector<long long> dist;
    vector<int> pre, pre_id;
    vector<bool> inq;
    void init(int n) { g.assign(n, {}); }
    void add_edge(int u, int v, long long cap, long long cost) {
        g[u].push_back({v, (int)g[v].size(), cap, cost});
        g[v].push_back({u, (int)g[u].size() - 1, 0, -cost});
    }
bool spfa(int s, int t) {
       dist.assign(g.size(), 1e18);
       inq.assign(g.size(), false);
       pre.assign(g.size(), -1);
       pre_id.assign(g.size(), -1);
       queue<int> q;
       dist[s] = 0; q.push(s); inq[s] = true;
       while (!q.empty()) {
           int u = q.front(); q.pop(); inq[u] = false;
           for (int i = 0; i < (int)g[u].size(); ++i) {
               auto &e = g[u][i];
               if (e.cap > 0 && dist[e.to] > dist[u] + e.cost) {
                   dist[e.to] = dist[u] + e.cost;
                   pre[e.to] = u;
                   pre_id[e.to] = i;
                   if (!inq[e.to]) {
                       q.push(e.to);
                       inq[e.to] = true;
                   }
               }
           }
       }
       return dist[t] != 1e18;
   }
pair<long long, long long> min_cost_flow(int s, int t) {
        long long flow = 0, cost = 0;
        while (spfa(s, t)) {
            long long f = 1e18;
            for (int v = t; v != s; v = pre[v])
                f = min(f, g[pre[v]][pre_id[v]].cap);
            flow += f;
            for (int v = t; v != s; v = pre[v]) {
                auto &e = g[pre[v]][pre_id[v]];
                e.cap -= f;
                g[v][e.rev].cap += f;
                cost += f * e.cost;
            }
        }
        return {flow, cost};
    }
};
如果费用非负且追求稳定复杂度，可用 Dijkstra + 势能 代替 SPFA。
```
= #text("开区间二分(l,r)")
```cpp
它的核心思想可以总结为“红蓝染色法”（维护不变量）。
💡 核心法则：我们将数组分为两半，一半是不满足条件的（红色，设为 l），一半是满足条件的（蓝色，设为 r）。
l 始终指向“不满足条件”的位置。r 始终指向“满足条件”的位置。
1. 初始化：因为 l 必须一开始就不满足条件，r 必须一开始就满足条件，我们要把它们放在数组界外（或者已知的边界上）。对于一个长度为 m 的数组 b（下标 0 - m - 1）
int l = -1; // -1 绝对不满足条件（越界了，算作不满足）
int r = m;  // m 绝对满足条件（我们要找的答案一定在 m 左侧，所以 m 兜底）
2. 循环条件：既然是开区间 (l, r)，说明 l 和 r 之间还有元素。只要它们俩之间还有至少一个元素，就继续二分
while (l + 1 < r) // 当 l 和 r 相邻时（比如 l=2, r=3），循环结束
3. 状态转移（极度舒适，不需要 +1 或 -1）
int mid = l + (r - l) / 2;
if (check(mid) == true) {
    r = mid; // mid 满足条件，所以把"蓝边界" r 移到 mid
} else {
    l = mid; // mid 不满足条件，所以把"红边界" l 移到 mid
}
4. 退出循环时的特点（重点！）：退出时，必定有 l + 1 == r（即 l 和 r 紧紧挨在一起）。因为 l 始终是不满足条件的最后一个位置，r 始终是满足条件的第一个位置。所以，你要找的目标答案就是 r！ 根本不用额外去推导什么边界。
```
= #text("线性求逆元")
```cpp
inv_num[1] = 1;
for (int i = 2; i < limit; i++) {
    inv_num[i] = (mod - mod / i) * inv_num[mod % i] % mod;
}
```
= #text("数位dp")
```cpp
string S;  // 上限
int n = S.size();
// dp[pos][state][tight]
// state 按题目定义（余数/和/积等）
vector<vector<vector<int>>> dp(n + 1, vector<vector<int>>(MAX_STATE, vector<int>(2, 0)));
dp[0][初始状态][1] = 1;  // 第0位，初始状态，受限
for (int i = 0; i < n; i++) {
    for (int st = 0; st < MAX_STATE; st++) {
        for (int tight = 0; tight < 2; tight++) {
            if (dp[i][st][tight] == 0) continue;

            int limit = tight ? (S[i] - '0') : 9;
            for (int d = 0; d <= limit; d++) {
                int nst = 转移(st, d);          // 按题目改
                int ntight = tight && (d == limit);
                dp[i + 1][nst][ntight] = (dp[i + 1][nst][ntight] + dp[i][st][tight]) % MOD;
            }
        }
    }
}
// 最终答案：所有 tight 状态的和
for (int tight = 0; tight < 2; tight++)
    ans = (ans + dp[n][目标状态][tight]) % MOD;
string S;  // 上限数字的字符串
int D;     // 题目参数
int n;
int memo[10005][105][2][2];  // 按题目需求调整大小
int dfs(int pos, int sum, bool tight, bool leadZero) {
    if (pos == n) {
        return sum == 0 ? 1 : 0;  // 终止条件按题目改
    }
    if (memo[pos][sum][tight][leadZero] != -1) return memo[pos][sum][tight][leadZero];
    int res = 0;
    int limit = tight ? (S[pos] - '0') : 9;
    for (int d = 0; d <= limit; d++) {
        bool ntight = tight && (d == limit);
        bool nleadZero = leadZero && (d == 0);

        if (leadZero && d == 0) {
            res = (res + dfs(pos + 1, sum, ntight, true)) % MOD;  // 前导零，不贡献
        } else {
            int nsum = (sum + d) % D;  // 状态更新按题目改
            res = (res + dfs(pos + 1, nsum, ntight, false)) % MOD;
        }
    }
    return memo[pos][sum][tight][leadZero] = res;
}
```
= #text("关于gcd的常数优化")
```cpp
int gcd(int a, int b) {
    int az = __builtin_ctz(a);
    int bz = __builtin_ctz(b);
    int z = min(az, bz);
    b >>= bz;
    while (a) {
        a >>= az;
        int diff = a - b;
        az = __builtin_ctz(diff);
        b = min(a, b), a = abs(diff);
    }
    return b << z;
}
```

```cpp
vector<array<int,3>> cal_(int x) {
    vector<array<int,3>> re(x + 1, {1, 1, 1});
    vector<int> not_prime(x + 1, 0);
    vector<int> prime;
    prime.reserve(x / 10);
    for (int i = 2; i <= x; i++) {
        if (!not_prime[i]) {
            prime.emplace_back(i);
            re[i] = {1, 1, i};
        }
        for (int j : prime) {
            if (i * j > x) break;
            not_prime[i * j] = 1;
            auto &[a, b, c] = re[i];
            int nj = a * j;
            if (nj <= b) {
                re[i * j] = {nj, b, c};
            } else if (nj <= c) {
                re[i * j] = {b, nj, c};
            } else {
                re[i * j] = {b, c, nj};
            }

            if (i % j == 0) break;
        }
    }
    return re;
}
int fast_gcd(int a, int b) {
    if (!a || !b) return (a | b);
    int az = __builtin_ctz(a);
    int bz = __builtin_ctz(b);
    int z = min(az, bz);
    b >>= bz;
    while(a) {
        a >>= az;
        int diff = a - b;
        az = __builtin_ctz(diff);
        b = min(a, b);
        a = abs(diff);
    }
    return b << z;
}
vector<vector<int>> get_gcd(int x) {
    int up = 1e3;
    vector<vector<int>> re(up + 1, vector<int>(up + 1));
    for (int i = 1; i <= up; i++) {
        for (int j = 1; j <= i; j++) {
            re[i][j] = re[j][i] = fast_gcd(i, j);
        }
    }
    for (int i = 0; i <= up; i++) {
        re[0][i] = re[i][0] = i;
    }
    return re;
}
/*
            auto &[x, y, z] = fact[b[j]];
            int _i = a[i];
            int gcd;
            if (x > 1) {
                gcd = _gcd[_i % x][x];
                tmp = tmp * gcd % mod;
                _i /= gcd;
            }

            if (y > 1) {
                gcd = _gcd[_i % y][y];
                tmp = tmp * gcd % mod;
                _i /= gcd;
            }

            if (z > 1e3) {
                tmp = tmp * (_i % z == 0 ? z : 1) % mod;
            }
            else {
                gcd = _gcd[_i % z][z];
                tmp = tmp * gcd % mod;
            }
*/
```
= #text("高精度")
#text("（1）版本一（长功能多）")

```cpp
struct BigInt {
    static const int BASE = 1e9;
    static const int WIDTH = 9;
    vector<int> s;

    BigInt(long long num = 0) { *this = num; }
    BigInt(string str) { *this = str; }

    BigInt& operator=(long long num) {
        s.clear();
        do { s.push_back(num % BASE); num /= BASE; } while (num > 0);
        return *this;
    }

    BigInt& operator=(const string& str) {
        s.clear();
        for (int i = str.length(); i > 0; i -= WIDTH) {
            if (i < WIDTH) s.push_back(stoi(str.substr(0, i)));
            else s.push_back(stoi(str.substr(i - WIDTH, WIDTH)));
        }
        trim();
        return *this;
    }

    void trim() {
        while (s.size() > 1 && s.back() == 0) s.pop_back();
    }

    BigInt operator+(const BigInt& b) const {
        BigInt c; c.s.clear();
        for (int i = 0, carry = 0; i < s.size() || i < b.s.size() || carry; ++i) {
            if (i < s.size()) carry += s[i];
            if (i < b.s.size()) carry += b.s[i];
            c.s.push_back(carry % BASE);
            carry /= BASE;
        }
        return c;
    }

    BigInt operator-(const BigInt& b) const {
        BigInt c = *this;
        for (int i = 0, borrow = 0; i < c.s.size(); ++i) {
            borrow = c.s[i] - borrow - (i < b.s.size() ? b.s[i] : 0);
            c.s[i] = borrow < 0 ? borrow + BASE : borrow;
            borrow = borrow < 0 ? 1 : 0;
        }
        c.trim();
        return c;
    }

    BigInt operator*(const BigInt& b) const {
        BigInt c; c.s.assign(s.size() + b.s.size(), 0);
        for (int i = 0; i < s.size(); ++i) {
            long long carry = 0;
            for (int j = 0; j < b.s.size() || carry; ++j) {
                long long cur = c.s[i + j] + s[i] * 1ll * (j < b.s.size() ? b.s[j] : 0) + carry;
                c.s[i + j] = cur % BASE;
                carry = cur / BASE;
            }
        }
        c.trim();
        return c;
    }

    bool operator<(const BigInt& b) const {
        if (s.size() != b.s.size()) return s.size() < b.s.size();
        for (int i = s.size() - 1; i >= 0; --i)
            if (s[i] != b.s[i]) return s[i] < b.s[i];
        return false;
    }

    bool operator==(const BigInt& b) const {
        return !(*this < b) && !(b < *this);
    }

    friend istream& operator>>(istream& in, BigInt& x) {
        string str; if (in >> str) x = str; return in;
    }

    friend ostream& operator<<(ostream& out, const BigInt& x) {
        out << x.s.back();
        for (int i = x.s.size() - 2; i >= 0; --i) {
            out << setfill('0') << setw(WIDTH) << x.s[i];
        }
        return out;
    } };
```

#text("（2）版本二（较短，功能简单）")

```cpp
struct Big {
    vector<int> v;
    const int B = 1e9;

    Big() {} // 默认空向量
    Big(long long n) { do { v.push_back(n % B); n /= B; } while (n); }
    Big(string s) {
        // 利用 max/min 优雅避开边界 if/else 截取判定
        for (int i = s.size(); i > 0; i -= 9)
            v.push_back(stoi(s.substr(max(0, i - 9), min(i, 9))));
        trim();
    }

    void trim() { while (v.size() > 1 && v.back() == 0) v.pop_back(); }

    Big operator+(const Big &b) const {
        Big c;
        // 进位 t 融合进循环条件，省去最后的 if(t) 判定
        for (int i = 0, t = 0; i < v.size() || i < b.v.size() || t; ++i) {
            if (i < v.size()) t += v[i];
            if (i < b.v.size()) t += b.v[i];
            c.v.push_back(t % B);
            t /= B;
        }
        return c;
    }

    Big operator-(const Big &b) const {
        Big c = *this;
        for (int i = 0, t = 0; i < c.v.size(); ++i) {
            t = c.v[i] - t - (i < b.v.size() ? b.v[i] : 0);
            c.v[i] = (t + B) % B;
            t = t < 0; // 核心：t < 0 时自动转为 1 形成借位，否则为 0
        }
        c.trim(); return c;
    }
Big operator*(const Big &b) const {
        Big c; c.v.assign(v.size() + b.v.size(), 0);
        for (int i = 0; i < v.size(); ++i) {
            long long t = 0;
            // 乘法的进位累加直接一气呵成
            for (int j = 0; j < b.v.size() || t; ++j) {
                t += c.v[i + j] + 1ll * v[i] * (j < b.v.size() ? b.v[j] : 0);
                c.v[i + j] = t % B;
                t /= B;
            }
        }
        c.trim(); return c;
    }

    bool operator<(const Big &b) const {
        if (v.size() != b.v.size()) return v.size() < b.v.size();
        for (int i = v.size() - 1; i >= 0; --i)
            if (v[i] != b.v[i]) return v[i] < b.v[i];
        return false;
    }

    friend ostream& operator<<(ostream& os, const Big& a) {
        if (a.v.empty()) return os << 0;
        os << a.v.back();
        for (int i = a.v.size() - 2; i >= 0; --i)
            os << setfill('0') << setw(9) << a.v[i];
        return os;
    }
};
```
= #text("闭区间线段树")
```cpp
    主要是文字说明一些注意点，比如左右区间为[l,m]和[m,r]，tl<m，走左，tr>m走右，点是没有长度的，长度由两个点相减得到，叶子节点r-l==1，r-l>=1才算有效区间，在一些区间问题常用这个写法，比如线段长度
```
#text("待施工： fwt，博弈论，根号分治，调和级数，点分治，polya定理带权重的版本，猫树，无旋treap，splay树，区间gcd最多下降log次，斐波那契数列的性质应用（每项大等于前一项，每一项小等于前一项的两倍，每一项等于前两项的和），重心点分治每次规模除二，st表二分，对顶堆，pbds，isap被特意卡的话常数比dinic大，")

#text("// 1. (x&b)>(b>>1) 判断x在b的最高位是否为1")

#text("// a + b == (a ^ b) + 2 * (a&b)")

#text("std::set：通过 Key、Compare、Allocator 三个参数，定义键类型、排序规则和内存分配。")

#text("std::map：在 set 的基础上增加 T 值类型，形成键、值、比较器、分配器四个可调参数。")

#text("std::unordered_set：用 Hash、KeyEqual 替代 Compare，组合成键、哈希、判等、分配器四个参数的无序集合。")

#text("std::unordered_map：在 unordered_set 模板参数上增加 T 值类型，成为键、值、哈希、判等、分配器五个可调参数的无序映射。")

#text("Bfs trick : Problem - B - Codeforces（来自cc cpc）")

#text("四个边框加入队列中bfs，可方便解决左上角和右下角是否连通等问题（本题左边框和下边框作为起点bfs，扫八个方向，标记为1，右边框与上边框作为起点，扫八个方向，标记为2，如果1与2连接（第二次八方扫描时扫描到），则左上与右下被隔绝）")

#text("当前树平均深度为Am，则加入一点作为原图上点的新儿子的期望深度为Am + 1,有时候期望计算可以从平均角度考虑")

#text("二进制分组")

#text("根据经典结论全零序列区间  得到序列  需要花费")

#text("数与数之间的关系问题转化为图上问题")

#text("矩阵乘法优化dp 问题")

#text("lcp优化--Ukkonen 算法或 Myers 差分算法")

#text("期望dp，计数dp（插入dp），动态dp，状压dp")

#text("Sub = (sub – 1) & mask // 遍历该集合的所有子集")

#text("笛卡尔树，差分约束，猫树，折半搜索，判断能否构成立体图形计算某点所有距离平方")

