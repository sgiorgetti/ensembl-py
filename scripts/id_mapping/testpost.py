import requests, sys
 
# server = "https://rest.ensembl.org"
# ext = "/lookup/id"
# headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
# r = requests.post(server+ext, headers=headers, data='{ "ids" : ["ENSG00000157764", "ENSG00000248378" ] }')

server = "http://registry-grpc:8080/registry/metaSearch"
headers={ "Content-Type" : "application/json" }
mydata={
  "name_pattern":"homo_sapiens_core_10%",
   "filters":[
      {
         "meta_key":"assembly.name",
         "meta_value":"GRCh38.p13"
      }
   ],
   "servers":[
      "anonymous@ensembldb.ensembl.org:3306"
   ]
}
r = requests.post(server, headers=headers, json=mydata)
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
print(repr(decoded.get('anonymous@ensembldb.ensembl.org:3306')[-1]))