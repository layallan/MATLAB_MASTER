import numpy as np
import matplotlib as plt

def loadData(flieName):

    lnum=0
    x=[]  
    
    with open(flieName, 'r') as f:
        for line in f:
            lnum += 1
            if(lnum>=16): 
                line=line.strip('\n') 
                            
                x.append(line[0])
                
 
    x=np.array(x)  
    # x=x.astype(np.float).tolist()  
    
    return x

def plotdata(x,y):
    fig = plt.figure(figsize=(10, 10))  
    ax = fig.add_subplot(1,1,1)
    ax.plot(x, y, 'red', label='unknown')  
    ax.legend(loc='upper left')  
    ax.set_xlabel('x-axis')  
    ax.set_ylabel('y-axis')  
   
    plt.show()  

if __name__ == '__main__':
    doc = r"5-12.txt"
    time = []
    voltage = []
    voltage = loadData(doc)
    print (voltage)