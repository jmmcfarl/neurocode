function MB = size_MB(X)

sz = whos('X');
MB = sz.bytes/1e6;