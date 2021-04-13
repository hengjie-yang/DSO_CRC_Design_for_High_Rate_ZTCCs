% A simple min-heap based PriorityQueue, which takes [a,b] as input and use
% b as the comparator to sort the element in the queue.
%
%
%Hengjie Yang - - Fall 2019
%
%
%Supported methods
% 
%       q = PriorityQueue();
%       q.insert([63,1]);
%       q.remove();
%       q.peek();
%       q.clear();

classdef PriorityQueue < handle
    properties (Access = private)
        Data % the data = [a, b] in the queue
        Size % the size of the queue
    end
    methods
        % Constructor 
        function obj = PriorityQueue(varargin)
            obj.Data = {};
            obj.Size = 0;
            if size(varargin) == 1
                [row, col] = size(varargin{1});
                if row == 1 && col == 2
                    obj.Data = varargin;
                    obj.Size = 1;
                end
            end
        end
        
            
        
        
        % insert a pair [a, b] into the queue given by arg1
        function insert(obj, arg1)
            if size(arg1) ~= 2
                disp('error: inserted element is invalid!');
                return
            end
            obj.Size = obj.Size + 1;
            obj.Data{obj.Size} = arg1;
            parent = floor(obj.Size/2);
            current = obj.Size;
            
            %perform bubble up operation
            while parent > 0
               if obj.Data{current}(2) < obj.Data{parent}(2)
                   % the child is smaller than its parent, swap!
                   temp = obj.Data{current};
                   obj.Data{current} = obj.Data{parent};
                   obj.Data{parent} = temp;
                   current = parent;
                   parent = floor(current/2);
               else
                  break;
               end
            end
        end
        
        
        
        %remove and return the top element in the queue
        function node = remove(obj)
            node = [Inf,-1];
            if obj.Size == 0
                disp('Remove from an empty queue, returning [Inf, -1]');
                return;
            end
            node = obj.Data{1};
            obj.Data{1} = obj.Data{obj.Size};
            obj.Data(obj.Size) = [];
            obj.Size = obj.Size - 1;
            current = 1;
            left = 2;
            right = 3;
            %perform bubble down
            while current < obj.Size
                if left <= obj.Size && obj.Data{left}(2) < obj.Data{current}(2)
                    smallest = left;
                else
                    smallest = current;
                end
                if right <= obj.Size && obj.Data{right}(2) < obj.Data{smallest}(2)
                    smallest = right;
                end
                if smallest ~= current
                    temp = obj.Data{current};
                    obj.Data{current} = obj.Data{smallest};
                    obj.Data{smallest} = temp;
                    current = smallest;
                    left = current*2;
                    right = current*2 + 1;
                else
                    break
                end
            end
            
        end
        
        
        %peek and return the first element of the queue
        function node = peek(obj)
            if obj.Size == 0
                node = [Inf, -1];
                disp('Queue is empty, returning [Inf, -1]');
                return
            end
            node = obj.Data{1};
        end
        
        
        %return the size of the queue
        function sz = size(obj)
            sz = obj.Size;
        end
        
        
        
        %clear the entire queue
        function clear(obj)
            obj.Data = {};
            obj.Size = 0;
        end
        
        
        %return all elements in the queue
        function queue = elements(obj)
            queue = obj.Data;
        end
            
           
    end
end
        
        

