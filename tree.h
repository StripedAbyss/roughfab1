#pragma once
#include<vector>
#include<iostream>
#include <fstream>

// 定义树节点类
class TreeNode {
public:
    int id; // 节点标识符
    std::vector<TreeNode*> children; // 子节点列表

    // 构造函数，初始化节点标识符
    TreeNode(int tid) : id(tid) {}
    void generateNodeGML(std::ofstream& outfile);
    void generateEdgeGML(std::ofstream& outfile);
    void Convert_GmlTree_To_GML();
};

// 定义树类
class Tree {
public:
    TreeNode* root; // 根节点指针

    // 构造函数，初始化根节点为空
    Tree() : root(nullptr) {}

    void insert(int parentId, int Id);
    void insertRec(TreeNode* node, int parentId, int Id);
    void printTree(TreeNode* node);
};