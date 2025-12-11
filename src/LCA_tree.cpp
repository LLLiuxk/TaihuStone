#include "modelGen.h"

void LCA_Tree::build(const AdjacencyList& adj)
{
    tree = adj;
    n = tree.size();
    LOG = std::ceil(std::log2(n));

    parent.assign(LOG + 1, std::vector<int>(n, -1));
    depth.assign(n, 0);

    // 以 0 为根进行 BFS
    int root = 0;

    std::queue<int> q;
    q.push(root);
    parent[0][root] = -1;
    depth[root] = 0;

    while (!q.empty()) {
        int u = q.front(); q.pop();

        for (int v : tree[u]) {
            if (v == parent[0][u]) continue;

            parent[0][v] = u;
            depth[v] = depth[u] + 1;
            q.push(v);
        }
    }

    // DP 填 parent[k][v]
    for (int k = 1; k <= LOG; k++) {
        for (int v = 0; v < n; v++) {
            int mid = parent[k - 1][v];
            parent[k][v] = (mid == -1 ? -1 : parent[k - 1][mid]);
        }
    }
}

int LCA_Tree::lca(int u, int v) const
{
    if (depth[u] < depth[v])
        std::swap(u, v);

    int diff = depth[u] - depth[v];

    // 提升 u
    for (int k = 0; k <= LOG; k++) {
        if (diff & (1 << k))
            u = parent[k][u];
    }

    if (u == v)
        return u;

    // 二分提升
    for (int k = LOG; k >= 0; k--) {
        if (parent[k][u] != parent[k][v]) {
            u = parent[k][u];
            v = parent[k][v];
        }
    }

    return parent[0][u];
}

std::vector<int> LCA_Tree::get_path(int u, int v) const
{
    int w = lca(u, v);

    std::vector<int> up, down;

    // u -> w
    int cur = u;
    while (cur != w) {
        up.push_back(cur);
        cur = parent[0][cur];
    }
    up.push_back(w);

    // v -> w
    cur = v;
    while (cur != w) {
        down.push_back(cur);
        cur = parent[0][cur];
    }
    // down = v -> ... -> w，需反转
    std::reverse(down.begin(), down.end());

    // 合并
    up.insert(up.end(), down.begin(), down.end());
    return up;
}