public void put(int key) {
    root = put(root, key);
    if(isRed(root)) {
        numRedNodes -=1;
    }
    root.color = BLACK;
}

// insert the key-value pair in the subtree rooted at h
private Node put(Node h, int key) {
    if (h == null) {
        numRedNodes++;
        return new Node(key, RED, 1);
    }

    int cmp = key - h.key;
    if      (cmp < 0) h.left  = put(h.left,  key);
    else if (cmp > 0) h.right = put(h.right, key);
    else              h.key   = key;

    if (!isRed(h.right) && isRed(h.left)) {
        h = rotateRight(h);
    }
    if (isRed(h.right)  &&  isRed(h.right.right)) {
        h = rotateLeft(h);
    }
    if (isRed(h.left)  &&  isRed(h.right)) {
        flipColors(h);
    }
    h.size = size(h.left) + size(h.right) + 1;

    return h;
}
