#include "modelGen.h"

PathQuery::PathQuery(int num_nodes, const std::vector<std::vector<int>>& adj, int root) 
    : n(num_nodes), adj_list(adj), root(root), parent(num_nodes, -1)
{
    isConnected = build_parent_tree();
    //if(root == 145) {
    //    for (int i = 0; i < num_nodes; i++)
    //        cout << "i: " << i << "  " << parent[i] << endl;
	//}
}

std::vector<int> PathQuery::query_path(int t) 
{
    if (t < 0 || t >= n) return {};
    if (parent[t] == -1 && t != root) {
        // 不连通，返回空路径
        return {};
    }
    std::vector<int> path;
    int cur = t;
    while (cur != -1) {
        path.push_back(cur);
        if (cur == root) break;
        cur = parent[cur];
    }
    std::reverse(path.begin(), path.end());
    return path;
}

// 一次 BFS 建立 parent 数组
bool PathQuery::build_parent_tree() 
{
    std::vector<bool> visited(n, false);
    std::queue<int> q;

    visited[root] = true;
    q.push(root);
    parent[root] = -1;
    int count = 1;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int v : adj_list[u]) {
            if (!visited[v]) {

                visited[v] = true;
                parent[v] = u;  // 记录路径树
                q.push(v);
                count++;
            }
        }
    }
    return count == n;
}