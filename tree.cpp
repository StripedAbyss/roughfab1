#include "tree.h"

// 生成GML文件中节点的部分
void TreeNode::generateNodeGML(std::ofstream& outfile) {
    // 写入节点的开始标记和ID到输出文件中
    outfile << "   node" << std::endl;
    outfile << "   [" << std::endl;
    outfile << "      id " << id << std::endl;
    outfile << "   ]" << std::endl;
    // 递归处理每个子节点
    for (TreeNode* child : children) {
        child->generateNodeGML(outfile);
    }
}

// 生成GML文件中边的部分
void TreeNode::generateEdgeGML(std::ofstream& outfile)
{
    // 遍历当前节点的每个子节点
    for (TreeNode* child : children) {
        // 写入边的开始标记以及源节点和目标节点到输出文件中
        outfile << "   edge" << std::endl;
        outfile << "   [" << std::endl;
        outfile << "      source " << id << std::endl;
        outfile << "      target " << child->id << std::endl;
        outfile << "   ]" << std::endl;
        // 递归处理当前子节点的子节点
        child->generateEdgeGML(outfile);
    }
}

// 将GmlTree转换为GML格式文件
void TreeNode::Convert_GmlTree_To_GML() {
    // 打开文件以写入GML数据
    std::ofstream outfile("tree.gml");
    if (outfile.is_open()) {
        // 写入GML文件的头部信息，包括图的类型
        outfile << "graph" << std::endl;
        outfile << "[" << std::endl;
        outfile << "   directed 0" << std::endl;

        // 生成GML文件中节点的部分
        generateNodeGML(outfile);
        // 生成GML文件中边的部分
        generateEdgeGML(outfile);

        // 写入GML文件的尾部信息，关闭文件
        outfile << "]" << std::endl;
        outfile.close();
    }
    else {
        // 如果文件打开失败，则输出错误消息
        std::cerr << "Unable to open file for writing." << std::endl;
    }
}

// 插入节点方法，根据父节点标识符和新节点标识符插入新节点
void Tree::insert(int parentId, int Id) {
    // 如果根节点为空，且父节点标识符为-1，则创建根节点
    if (root == nullptr) {
        if (parentId == -1) {
            root = new TreeNode(Id);
        }
        return;
    }

    // 递归插入新节点
    insertRec(root, parentId, Id);
}

// 递归插入节点方法
void Tree::insertRec(TreeNode* node, int parentId, int Id) {
    // 如果当前节点的标识符等于父节点标识符，则将新节点添加为当前节点的子节点
    if (node->id == parentId) {
        node->children.push_back(new TreeNode(Id));
    }
    else {
        // 否则递归查找子节点并插入新节点
        for (TreeNode* child : node->children) {
            insertRec(child, parentId, Id);
        }
    }
}

// 打印树方法，从指定节点开始递归打印树的结构
void Tree::printTree(TreeNode* node) {
    if (node == nullptr) {
        return;
    }

    // 打印当前节点的标识符以及其所有子节点的标识符
    std::cout << node->id << ": ";
    for (TreeNode* child : node->children) {
        std::cout << child->id << " ";
    }
    std::cout << std::endl;

    // 递归打印每个子节点的树结构
    for (TreeNode* child : node->children) {
        printTree(child);
    }
}

