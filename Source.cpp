// Graph implementation a graph using array
// Graph is stored in matrix of size n*n
// Graph operations are performed using nodes, edges -- the adjecency matrix.
// Graph is weighted/unweighted
// Graph is directed/undirected
// Matrix weights are non-negative
// Functions interfaces are {isEmpty, displayGraph, addEdge, addDiEdge}

#include <iostream>
#include <stdexcept>
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////   STACK   /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Stack {
    int top;
    const int size;
    int* array;
};

Stack* createStack(const int arg_size)
{
    Stack* stack = new Stack{ -1, arg_size };
    stack->array = new int[arg_size];
    return stack;
}

// FUNCTION:    Is Empty.
// REQUIRES:    NOTHING.
// MODIFIES:    NOTHING.
// Executes:
//              If the {top} is equal to sentinel value:
//                  Stack is empty.
//              If passed all the above guards:
//                  Stack is not empty.
bool isEmpty(const Stack*);

// FUNCTION:    Is Full.
// REQUIRES:    NOTHING.
// MODIFIES:    NOTHING.
// EXECUTES:    
//              If the {top} is last possible element in the stack:
//                  Stack is full.
//              If passed all the above guards:
//                  Stack is not full.
bool isFull(const Stack*);

// FUNCTION:    Push.
// REQUIRES:    That the stack is not full.
// MODIFIES:    {Top} of the stack. The stack.
// EXECUTES:    
//              If the stack is full:
//                  Throw overflow exception.
//              If passed all the above guards:
//                  Add new element after the {top}
//                  Increment {top}
void push(const int, Stack*);

// FUNCTION:    Peek.
// REQUIRES:    That the stack is not empty.
// MODIFIES:    NOTHING.
// EXECUTES:
//              If the stack is empty:
//                  Throw underflow exception.
//              If passed all the above guards:
//                  Return value at the {top} of the queue.
int peek(const Stack*);

// FUNCTION:    Search.
// REQUIRES:    Nothing.
// MODIFIES:    Nothing.
// EXECUTES:
//              Compare each value -- in the stack -- from start till {top} -- with searched value:
//                  If any pair matches:
//                      The searched value exists in the stack.
//              If passed all the above guards:
//                  The searched value do not exists in the stack.
bool search(const int, const Stack*);

// FUNCTION:    Pop
// REQUIRES:    That the stack is not empty.
// MODIFIES:    {Top} of the stack. The stack.
// EXECUTES:
//              If the stack is empty:
//                  Throw underflow exception.
//              If passed all the above guards:
//                  Remove the {top} from the stack.
//                  Return the removed value.
int pop(Stack*);

bool isFull(const Stack* ptr_Stack) {
    return (ptr_Stack->top == (ptr_Stack->size - 1));
}

bool isEmpty(const Stack* ptr_Stack) {
    return (ptr_Stack->top == -1);
}

void push(const int arg_value, Stack* ptr_Stack) {
    if (isFull(ptr_Stack)) {
        throw overflow_error("Stack is already full. Cannot push()");
    }
    ptr_Stack->array[++ptr_Stack->top] = arg_value;
    return;
}

int peek(const Stack* ptr_Stack) {
    if (isEmpty(ptr_Stack)) {
        throw underflow_error("Stack is already empty. Cannot peek()");
    }
    return ptr_Stack->array[ptr_Stack->top];
}

bool search(const int arg_value, const Stack* ptr_Stack) {
    for (int i = 0; i <= ptr_Stack->top; ++i) {
        if (arg_value == ptr_Stack->array[i]) {
            return(true);
        }
    }
    return(false);
}

int pop(Stack* ptr_Stack) {
    if (isEmpty(ptr_Stack)) {
        throw underflow_error("Stack is already empty. Cannot pop()");
    }
    return ptr_Stack->array[ptr_Stack->top--];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////   QUEUE   /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Queue {
    int head;
    int tail;
    const int size;
    int* array;
    int& front = head;
    int& rear = tail;
};

Queue* createQueue(const int arg_size)
{
    Queue* queue = new Queue{ -1, -1, arg_size };
    queue->array = new int[arg_size];
    return queue;
}

// FUNCTION:    Is Empty.
// REQUIRES:    NOTHING.
// MODIFIES:    NOTHING.
// EXECUTES:    
//              If the {head/front} is equal to {sentinel value}:
//                  Queue is empty.
//              If passed all the above guards:
//                  Queue is not empty.
bool isEmpty(const Queue* ptr_Queue);

// FUNCTION:    Is Full.
// REQUIRES:    NOTHING.
// MODIFIES:    NOTHING.
//              Calculate the {next position} -- for queue -- of insertion.
//              If the {next position} is the {head/front} of the queue:
//                  Queue is full.
//              If passed all the above guards:
//                  Queue is not full.
bool isFull(const Queue* ptr_Queue);

// FUNCTION:    Enqueue.
// REQUIRES:    That the queue is not full.
// MODIFIES:    Rear, and front of the queue. The queue.
// EXECUTES:    
//              If the queue is full:
//                  Throw overflow exception.
//              If the queue is empty:
//                  Insert the new element as first element of the queue.
//                  This new element will be the {head/front} as well as {tail/rear}
//              If passed all the above guards:
//                  Add new element after the {tail/rear}
//                  Increment {tail/rear}
void enqueue(const int, Queue* ptr_Queue);

// FUNCTION:    Peek.
// REQUIRES:    That the queue is not empty.
// MODIFIES:    NOTHING.
// EXECUTES:
//              If the queue is empty:
//                  Throw underflow exception.
//              If passed all the above guards:
//                  Return value at the {head/front} of the queue.
int peek(const Queue* ptr_Queue);

// FUNCTION:    Search.
// REQUIRES:    NOTHING.
// MODIFIES:    NOTHING.
// EXECUTES:
//              Compare each value starting -- in the queue -- from the {head/front} to {tail/rear} -- with searched value:
//                  If any pair matches:
//                      The searched value exists in the queue.
//              If passed all the above guards:
//                  The searched value do not exists in the queue.
bool search(const int, const Queue* ptr_Queue);

// FUNCTION:    Dequeue
// REQUIRES:    That the queue is not empty.
// MODIFIES:    Rear, and front of the queue. The queue.
// EXECUTES:
//              If the queue is empty:
//                  Throw underflow exception.
//              If their is only one element in the queue:
//                  Remove that element from the queue.
//                  Reset the {head/front} and {tail/rear} to sentinel value.
//                  Return the removed value.
//              If passed all the above guards:
//                  Remove the {head/front} from the queue.
//                  Return the removed value.
int dequeue(Queue* ptr_Queue);

bool isEmpty(const Queue* ptr_Queue) {
    return(ptr_Queue->front == -1);
}

bool isFull(const Queue* ptr_Queue) {
    const int nextPos = ((ptr_Queue->rear + 1) % ptr_Queue->size);
    return(nextPos == ptr_Queue->front);
}

int getElements(const Queue* ptr_Queue) {
    if (isEmpty(ptr_Queue))
        return 0;
    if (ptr_Queue->rear >= ptr_Queue->front) {
        return ((ptr_Queue->rear - ptr_Queue->front) + 1);
    }
    return ((ptr_Queue->rear + 1) + (ptr_Queue->size - ptr_Queue->front));
}

void enqueue(const int arg_value, Queue* ptr_Queue) {
    if (isFull(ptr_Queue)) {
        throw overflow_error("Queue is already full. Cannot enqueue()");
    }
    if (isEmpty(ptr_Queue)) {
        ++(ptr_Queue->front);
        ptr_Queue->array[++(ptr_Queue->rear)] = arg_value;
        return;
    }
    ptr_Queue->array[++(ptr_Queue->rear)] = arg_value;
    return;
}

int peek(const Queue* ptr_Queue) {
    if (isEmpty(ptr_Queue)) {
        throw overflow_error("Queue is already empty. Cannot peek()");
    }
    return ptr_Queue->array[ptr_Queue->front];
}

bool search(const int arg_value, const Queue* ptr_Queue) {
    for (int i = ptr_Queue->front; i <= ptr_Queue->rear; ++i) {
        if (arg_value == ptr_Queue->array[i]) {
            return(true);
        }
    }
    return(false);
}

int dequeue(Queue* ptr_Queue) {
    if (isEmpty(ptr_Queue)) {
        throw overflow_error("Queue is already empty. Cannot dequeue()");
    }
    if (ptr_Queue->rear == ptr_Queue->front)
    {
        int toReturn = ptr_Queue->array[ptr_Queue->front];
        ptr_Queue->rear = -1;
        ptr_Queue->front = -1;
        return toReturn;
    }
    return ptr_Queue->array[(ptr_Queue->front)++];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////   GRAPH   /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Graph {
    const int numberOfNodes;
    const bool isDirected;
    const bool isWeighted;
    int* nodes;
    int** adjMatrix;
};

Graph* createGraph(const int arg_numberOfNodes, const bool arg_isDirected, const bool arg_isWeighted) {
    Graph* graph = new Graph{ arg_numberOfNodes, arg_isDirected, arg_isWeighted };
    graph->nodes = new int[arg_numberOfNodes];
    for (int i = 0; i < (arg_numberOfNodes); ++i) {
        graph->nodes[i] = (i + 1);
    }
    graph->adjMatrix = new int* [arg_numberOfNodes];
    for (int i = 0; i < (arg_numberOfNodes); ++i) {
        graph->adjMatrix[i] = new int[arg_numberOfNodes];
        for (int j = 0; j < (arg_numberOfNodes); ++j) {
            graph->adjMatrix[i][j] = 0;
        }
    }
    return graph;
}

bool isEmpty(const Graph*);
int getGraphNodes(const Graph*);
int getGraphEdges(const Graph*);

int getGraphSinkNodes(const Graph*);
int getGraphSourceNodes(const Graph*);
int getGraphIsolatedNodes(const Graph*);
bool isNodeSink(const int, const Graph*);
bool isNodeSource(const int, const Graph*);
bool isNodeIsolated(const int, const Graph*);

int getNodeInDegree(const int, const Graph*);
int getNodeOutDegree(const int, const Graph*);
int getNodeTotalDegree(const int, const Graph*);

int getNodeNeighbourCount(const int, const Graph*);
Queue* getNodeNeighbours(const int, const Graph*);

Queue* getBreadthFirstTraversalNodes(const int, const Graph*);
bool areNodesConnected(const int, const int, const Graph*);
Queue* getSCCNodes(const int, const Graph*);
Queue* getWCCNodes(const int, const Graph*);
int getGraphMaxSCCCount(const Graph*);
int getGraphMaxWCCCount(const Graph*);
// Display Functions
void displayGraphEquation(const Graph*);
void displayGraphOutDegrees(const Graph*);
void displayGraphInDegrees(const Graph*);
void displayGraphTotalDegrees(const Graph*);
void displayGraphSinkNodes(const Graph*);
void displayGraphSourceNodes(const Graph*);
void displayGraphIsolatedNodes(const Graph*);
void displayGraphMaxSCCCount(const Graph*);
void displayGraphMaxWCCCount(const Graph*);
void displayGraphInfo(const Graph*);
void displayGraphMatrix(const Graph*);
void displaySCCDistribution(const Graph*);
void displayWCCDistribution(const Graph*);

void addEdge(const int, const int, Graph*);
void addEdge(const int, const int, const int, Graph*);

bool isEmpty(const Graph* ptr_graph) {
    if (ptr_graph->numberOfNodes == 0) {
        return true;
    }
    return false;
}

int getGraphNodes(const Graph* ptr_graph) {
    return ptr_graph->numberOfNodes;
}

int getGraphEdges(const Graph* ptr_graph) {
    int numberOfEdges = 0;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++) {
        for (int j = 0; j < ptr_graph->numberOfNodes; j++) {
            numberOfEdges += ptr_graph->adjMatrix[i][j];
        }
    }
    return numberOfEdges;
}

int getGraphSinkNodes(const Graph* ptr_graph) {
    int numberOfSinkNodes = 0;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++) {
        if (isNodeSink(i + 1, ptr_graph)) {
            numberOfSinkNodes += 1;
        }
    }
    return numberOfSinkNodes;
}

int getGraphSourceNodes(const Graph* ptr_graph) {
    int numberOfSourceNodes = 0;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++) {
        if (isNodeSource(i + 1, ptr_graph)) {
            numberOfSourceNodes += 1;
        }
    }
    return numberOfSourceNodes;
}

int getGraphIsolatedNodes(const Graph* ptr_graph) {
    int numberOfIsolatedNodes = 0;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++) {
        if (isNodeIsolated(i + 1, ptr_graph)) {
            numberOfIsolatedNodes += 1;
        }
    }
    return numberOfIsolatedNodes;
}

bool isNodeSink(const int arg_node, const Graph* ptr_graph) {
    if (arg_node < 1 || arg_node > getGraphNodes(ptr_graph)) {
        throw invalid_argument("Bad Nodes");
    }
    if (getNodeOutDegree(arg_node, ptr_graph) == 0 && getNodeInDegree(arg_node, ptr_graph) > 0) {
        return true;
    }
    else {
        return false;
    }

}

bool isNodeSource(const int arg_node, const Graph* ptr_graph) {
    if (arg_node < 1 || arg_node > getGraphNodes(ptr_graph)) {
        throw invalid_argument("Bad Nodes");
    }
    if (getNodeOutDegree(arg_node, ptr_graph) > 0 && getNodeInDegree(arg_node, ptr_graph) == 0) {
        return true;
    }
    else {
        return false;
    }

}

bool isNodeIsolated(const int arg_node, const Graph* ptr_graph) {
    if (arg_node < 1 || arg_node > getGraphNodes(ptr_graph)) {
        throw invalid_argument("Bad Nodes");
    }
    if (getNodeOutDegree(arg_node, ptr_graph) == 0 && getNodeInDegree(arg_node, ptr_graph) == 0) {
        return true;
    }
    else {
        return false;
    }

}

int getNodeOutDegree(const int arg_node, const Graph* ptr_graph) {
    if (arg_node < 1 || arg_node > getGraphNodes(ptr_graph)) {
        throw invalid_argument("Bad Nodes");
    }
    int i = arg_node - 1;
    int outDegree = 0;
    for (int j = 0; j < getGraphNodes(ptr_graph); j++) {
        if (ptr_graph->adjMatrix[i][j] != 0) {
            outDegree += ptr_graph->adjMatrix[i][j];
        }
    }
    return outDegree;
}

int getNodeInDegree(const int arg_node, const Graph* ptr_graph) {
    if (arg_node < 1 || arg_node > getGraphNodes(ptr_graph)) {
        throw invalid_argument("Bad Nodes");
    }
    int j = arg_node - 1;
    int outDegree = 0;
    for (int i = 0; i < getGraphNodes(ptr_graph); i++) {
        if (ptr_graph->adjMatrix[i][j] != 0) {
            outDegree += ptr_graph->adjMatrix[i][j];
        }
    }
    return outDegree;
}

int getNodeTotalDegree(const int arg_node, const Graph* ptr_graph) {
    int totalDegree = getNodeInDegree(arg_node, ptr_graph) + getNodeOutDegree(arg_node, ptr_graph);
    return totalDegree;
}

int getNodeNeighbourCount(const int arg_node, const Graph* ptr_graph) {
    if (arg_node < 1 || arg_node > getGraphNodes(ptr_graph)) {
        throw invalid_argument("Bad Nodes");
    }
    int i = arg_node - 1;
    int outDegree = 0;
    for (int j = 0; j < getGraphNodes(ptr_graph); j++) {
        if (ptr_graph->adjMatrix[i][j] != 0) {
            outDegree += 1;
        }
    }
    return outDegree;
}

Queue* getNodeNeighbours(const int arg_node, const Graph* ptr_graph) {
    if (arg_node < 1 || arg_node > getGraphNodes(ptr_graph)) {
        throw invalid_argument("Bad Nodes");
    }
    const int i = arg_node - 1;
    Queue* queue = createQueue(getNodeNeighbourCount(arg_node, ptr_graph));
    for (int j = 0; j < getGraphNodes(ptr_graph); ++j) {
        if (ptr_graph->adjMatrix[i][j] != 0)
        {
            enqueue(j + 1, queue);
        }
    }
    return queue;
}

Queue* getBreadthFirstTraversalNodes(const int arg_node, const Graph* ptr_graph) {
    if (arg_node < 1 || arg_node > getGraphNodes(ptr_graph)) {
        throw invalid_argument("Bad Nodes");
    }
    Queue* visited = createQueue(getGraphNodes(ptr_graph));
    Queue* mQueue = createQueue(getGraphNodes(ptr_graph));
    enqueue(arg_node, mQueue);
    while (!isEmpty(mQueue)) {
        Queue* neighbours = getNodeNeighbours(peek(mQueue), ptr_graph);
        while (!isEmpty(neighbours)) {
            if (search(peek(neighbours), mQueue) || search(peek(neighbours), visited)) {
                dequeue(neighbours);
            }
            else {
                enqueue(dequeue(neighbours), mQueue);
            }
        }
        if (search(peek(mQueue), visited)) {
            dequeue(mQueue);
        }
        else {
            enqueue(dequeue(mQueue), visited);
        }
        if (isFull(visited)) {
            return visited;
        }
    }
    return visited;
}

bool areNodesConnected(const int arg_fromNode, const int arg_toNode, const Graph* ptr_graph) {
    if (arg_fromNode < 1 || arg_fromNode > getGraphNodes(ptr_graph)) {
        throw invalid_argument("Bad Nodes");
    }
    if (arg_toNode < 1 || arg_toNode > getGraphNodes(ptr_graph)) {
        throw invalid_argument("Bad Nodes");
    }
    Queue* reachables = getBreadthFirstTraversalNodes(arg_fromNode, ptr_graph);
    while (!isEmpty(reachables)) {
        if (dequeue(reachables) == arg_toNode) {
            return true;
        }
    }
    return false;
}

Queue* getSCCNodes(const int arg_node, const Graph* ptr_graph) {
    if (arg_node < 1 || arg_node > getGraphNodes(ptr_graph)) {
        throw invalid_argument("Bad Nodes");
    }
    Queue* reachables = getBreadthFirstTraversalNodes(arg_node, ptr_graph);
    Queue* SCC = createQueue(getElements(reachables));
    while (!isEmpty(reachables)) {
        if (areNodesConnected(peek(reachables), arg_node, ptr_graph)) {
            enqueue(peek(reachables), SCC);
        }
        dequeue(reachables);
    }
    return SCC;
}

Queue* getWCCNodes(const int arg_node, const Graph* ptr_graph) {
    if (arg_node < 1 || arg_node > getGraphNodes(ptr_graph)) {
        throw invalid_argument("Bad Nodes");
    }
    Queue* reachables = getBreadthFirstTraversalNodes(arg_node, ptr_graph);
    Queue* WCC = createQueue(getElements(reachables));
    while (!isEmpty(reachables)) {
        if (!areNodesConnected(peek(reachables), arg_node, ptr_graph)) {
            enqueue(peek(reachables), WCC);
        }
        dequeue(reachables);
    }
    return WCC;
}

int getGraphMaxSCCCount(const Graph* ptr_graph) {
    int maxSCCCount = 1;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++) {
        Queue* queue = getSCCNodes(i + 1, ptr_graph);
        int count = getElements(queue);
        if (count > maxSCCCount) {
            maxSCCCount = count;
        }
    }
    return maxSCCCount;
}

int getGraphMaxWCCCount(const Graph* ptr_graph) {
    int maxWCCCount = 1;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++) {
        Queue* queue = getWCCNodes(i + 1, ptr_graph);
        int count = getElements(queue);
        if (count > maxWCCCount) {
            maxWCCCount = count;
        }
    }
    return maxWCCCount;
}

void displayGraphEquation(const Graph* ptr_graph) {
    cout << "G(" << getGraphNodes(ptr_graph) << ", " << getGraphEdges(ptr_graph) << ")" << endl;
}

void displayGraphOutDegrees(const Graph* ptr_graph) {
    cout << endl;
    cout << "Node" << "\tOut Degree" << endl;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++) {
        cout << ptr_graph->nodes[i] << "\t" << getNodeOutDegree(i + 1, ptr_graph) << endl;
    }
}
void displayGraphInDegrees(const Graph* ptr_graph) {
    cout << endl;
    cout << "Node" << "\tIn Degree" << endl;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++) {
        cout << ptr_graph->nodes[i] << "\t" << getNodeInDegree(i + 1, ptr_graph) << endl;
    }
}
void displayGraphTotalDegrees(const Graph* ptr_graph) {
    cout << endl;
    cout << "Node" << "\t\tIn" << "\t\tOut" << "\t\tTotal" << endl;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++) {
        cout << ptr_graph->nodes[i];
        cout << "\t\t" << getNodeOutDegree(i + 1, ptr_graph);
        cout << "\t\t" << getNodeInDegree(i + 1, ptr_graph);
        cout << "\t\t" << getNodeTotalDegree(i + 1, ptr_graph) << endl;
    }
}

void displayGraphSinkNodes(const Graph* ptr_graph) {
    cout << "Graph G has " << getGraphSinkNodes(ptr_graph) << " Sink Nodes" << endl;
}

void displayGraphSourceNodes(const Graph* ptr_graph) {
    cout << "Graph G has " << getGraphSourceNodes(ptr_graph) << " Source Nodes" << endl;
}

void displayGraphIsolatedNodes(const Graph* ptr_graph) {
    cout << "Graph G has " << getGraphIsolatedNodes(ptr_graph) << " Isolated Nodes" << endl;
}

void displayGraphMaxSCCCount(const Graph* ptr_graph) {
    cout << "Graph G has " << getGraphMaxSCCCount(ptr_graph) << " as Maximum SSC component" << endl;
}

void displayGraphMaxWCCCount(const Graph* ptr_graph) {
    cout << "Graph G has " << getGraphMaxWCCCount(ptr_graph) << " as Maximum SSC component" << endl;
}

void displayGraphInfo(const Graph* ptr_graph) {
    cout << "Calculating" << " " << "..." << endl << "\t\t";
    displayGraphEquation(ptr_graph);

    cout << "Calculating" << " GraphEquation " << "..." << endl << "\t\t";;
    displayGraphSinkNodes(ptr_graph);

    cout << "Calculating" << " SourceNodes " << "..." << endl << "\t\t";;
    displayGraphSourceNodes(ptr_graph);

    cout << "Calculating" << " IsolatedNodes " << "..." << endl << "\t\t";;
    displayGraphIsolatedNodes(ptr_graph);

    cout << "Calculating" << " SCCCount " << "..." << endl << "\t\t";;
    displayGraphMaxSCCCount(ptr_graph);

    cout << "Calculating" << " WCCCount " << "..." << endl << "\t\t";;
    displayGraphMaxWCCCount(ptr_graph);
}

void displayGraphMatrix(const Graph* ptr_graph) {
    cout << endl << endl;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++) {
        cout << '\t' << ptr_graph->nodes[i];
    }
    cout << endl;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++)
    {
        cout << ptr_graph->nodes[i];
        for (int j = 0; j < ptr_graph->numberOfNodes; j++)
        {
            cout << '\t' << ptr_graph->adjMatrix[i][j];
        }
        cout << endl;
    }
    return;
}

void displaySCCDistribution(const Graph* ptr_graph) {
    cout << endl << endl;
    cout << "Node" << "\t\tSCC" << "\t\tSCCs" << endl;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++) {
        cout << ptr_graph->nodes[i];
        Queue* queue = getSCCNodes(i + 1, ptr_graph);
        cout << "\t\t" << getElements(queue);
        cout << "\t\t";
        while (!isEmpty(queue)) {
            cout << dequeue(queue) << " ";
        }
        cout << endl;
    }
}

void displayWCCDistribution(const Graph* ptr_graph) {
    cout << endl << endl;
    cout << "Node" << "\t\tWCC" << "\t\tWCCs" << endl;
    for (int i = 0; i < ptr_graph->numberOfNodes; i++) {
        cout << ptr_graph->nodes[i];
        Queue* queue = getWCCNodes(i + 1, ptr_graph);
        cout << "\t\t" << getElements(queue);
        cout << "\t\t";
        while (!isEmpty(queue)) {
            cout << dequeue(queue) << " ";
        }
        cout << endl;
    }
}

void addEdge(const int arg_fromNode, const int arg_toNode, Graph* ptr_graph) {
    if (ptr_graph->isWeighted) {
        throw invalid_argument("Argument Missing: Weights. Passed no weights for an weighted graph.");
    }
    const int i = arg_fromNode - 1;
    const int j = arg_toNode - 1;
    if (ptr_graph->isDirected) {

        ptr_graph->adjMatrix[i][j] += 1;
        if (i == j)
            ptr_graph->adjMatrix[j][i] += 1;
        return;
    }
    else {
        ptr_graph->adjMatrix[i][j] += 1;
        ptr_graph->adjMatrix[j][i] += 1;
        return;
    }
}

void addEdge(const int arg_fromNode, const int arg_toNode, const int weight, Graph* ptr_graph) {
    if (!ptr_graph->isWeighted) {
        addEdge(arg_fromNode, arg_toNode, ptr_graph);
        return;
    }
    const int i = arg_fromNode - 1;
    const int j = arg_toNode - 1;
    if (ptr_graph->isDirected) {
        ptr_graph->adjMatrix[i][j] = weight;
        return;
    }
    else {
        ptr_graph->adjMatrix[i][j] = weight;
        ptr_graph->adjMatrix[j][i] = weight;
        return;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////   MAIN   //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <fstream>
string path = "./CA-GrQc.txt";
int main(void)
{
    ifstream file;
    file.open(path);
    if (!file.is_open()){
        cout << "File Not Opened. Exiting" << endl;
        exit(1);
    }
    int nodes;
    file >> nodes;
    Graph* mGraph = createGraph(nodes, true, false);

    int from;
    int to;
    while (file >> from)
    {
        file >> to;

        if ((from < 1 || from > nodes) || (to < 1 || to > nodes))
            continue;
        addEdge(from, to, mGraph);
    }

    displayGraphInfo(mGraph);
    //displayGraphMatrix(mGraph);
    //displayGraphTotalDegrees(mGraph);
    return 0;
}

//int main(void)
//{
//    Graph* mGraph = createGraph(8, true, false);
//
//    addEdge(1, 2, mGraph);
//    addEdge(2, 3, mGraph);
//    addEdge(3, 4, mGraph);
//    addEdge(3, 5, mGraph);
//    addEdge(4, 1, mGraph);
//    addEdge(5, 6, mGraph);
//    addEdge(6, 7, mGraph);
//    addEdge(7, 5, mGraph);
//    addEdge(7, 8, mGraph);
//
//    displayGraphInfo(mGraph);
//    displayGraphMatrix(mGraph);
//    displayGraphTotalDegrees(mGraph);
//    displaySCCDistribution(mGraph);
//    displayWCCDistribution(mGraph);
//
//    return 0;
//}