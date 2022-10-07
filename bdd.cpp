#include <iostream>
#include <map>
#include <vector>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

class BinaryDecisionDiagram
{
    public:
    struct Node
    {
        float lb, ub;
        Node *son0, *son1;
    };

    Node leaf0;
    Node leaf1;

    std::vector<std::vector<Node*>> tree;

    Node* buildQOBDD(int slack, int level, std::vector<int> &a, int b)
    {
        int s = 0;
        for(int i = level; i < a.size(); i++)
            s += a[i];

        if ((slack < 0) | (slack > s))
            return &leaf0;     
        
        if ((slack == 0) & (slack == s))
            return &leaf1;

        Node* u = new Node;
        tree[level].push_back(u);

        u->son0 = buildQOBDD(slack, level+1, a, b);
        u->son1 = buildQOBDD(slack - a[level], level+1, a, b);
        u->lb = max(u->son0->lb, u->son1->lb + a[level]);
        u->ub = min(u->son0->ub,u->son1->ub + a[level]);
        return u;
    }
    BinaryDecisionDiagram()
    {
        leaf0.lb = -999999;
        leaf0.ub = -1;
        leaf1.lb = 0;
        leaf0.ub = +999999;
        tree.resize(3);
    }
};

int main(int argc, char const *argv[])
{
    BinaryDecisionDiagram bdd = BinaryDecisionDiagram();
    std::vector<int> a = {1,1,1};
    int b = 1;
    bdd.buildQOBDD(b,0,a,b);

    for (int l = 0; l < 3; l++)
    {
        std::cout << "livello: " << l << std::endl;
        for (auto& n : bdd.tree[l])
        {
            std::cout << n << "\t[" << n->lb << "," << n->ub << "]\t" 
            << n->son0 << " " << n->son1 
            << std::endl;
        }
    }
    std::cout << "leaf0:" << &(bdd.leaf0) << "\n"
                << "leaf1:" << &(bdd.leaf1) << "\n";

    return 0;
}
