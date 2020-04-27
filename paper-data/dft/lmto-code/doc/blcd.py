



from pylab import *

f = figure(figsize=(20,10))




#def mkgrid(ax,p,d,ud):
   ##p = array(p)
   ##d = array(d)
   #for i in range(d[0]+1):
      #plot([p[0]+i*ud[0], p[0]+i*ud[0]], [p[1], p[1]+ud[1]*d[1]], ls='-', marker='')

   #for i in range(d[1]+1):
      #plot([p[0], p[0]+ud[0]*d[0]], [p[1]+i*ud[1], p[1]+i*ud[1]], ls='-', marker='')

#mkgrid(ax,[1,1], [16,16], [1,1])



rszs = array([25,25])
pszs = array([4,4])
bszs = array([3,3])
blks = (rszs - 1) // bszs + 1
lblk = (blks - 1) // pszs + 1
lszs = lblk * bszs
gszs = lszs * pszs

pcols = cm.get_cmap('Paired',pszs.prod())

ax = f.add_axes([0,0,1,1], aspect='equal', frameon=False)
ax.set_xticks([])
ax.set_yticks([])

ipos = array([1,1])
for i in range(blks[0]):
   for j in range(blks[1]):
      bcrd = array([i,j])
      pcrd = bcrd%pszs
      pid = pcrd[0]*pszs[1]+pcrd[1]
      llc = ipos + bcrd*bszs
      ax.add_patch(Rectangle(llc[::-1], bszs[1], bszs[0], fc=pcols(pid),ec='black'))
      ax.text((llc+bszs*0.5)[1],(llc+bszs*0.5)[0], \
         '%i,%i\n%i[%i,%i]'%(tuple(bcrd)+(pid,)+tuple(pcrd)), \
         ha='center', va='center')


ipos = array([1,ipos[1]+gszs[1]+2])
axp = ax
for i in range(pszs[0]):
   for j in range(pszs[1]):
      pcrd = array([i,j])
      llc = ipos + pcrd*lszs
      pid = pcrd[0]*pszs[1]+pcrd[1]
      ax.add_patch(Rectangle(llc[::-1], lszs[1], lszs[0], fc=pcols(pid),ec='black'))

      for mi in range(lblk[0]):
         for mj in range(lblk[1]):
            lcrd = array([mi,mj])

            mllc = llc + lcrd*bszs
            ax.add_patch(Rectangle(mllc[::-1], bszs[1], bszs[0], fc=pcols(pid), ec='black'))
            gcrd = lcrd*pszs+pcrd
            ax.text((mllc+0.5*bszs)[1],(mllc+0.5*bszs)[0], \
               '%i,%i'%(tuple(gcrd)), \
               ha='center', va='center')


ax.axis([0,2*(gszs[0]+2)+5,0,gszs[1]+2])
ax.invert_yaxis()

#show()
savefig('blockcyclic2D.pdf')