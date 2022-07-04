module github.com/tidwall/geodesic

go 1.17

retract (
	v1.52.3 // equal to v0.3.1
	v1.52.2 // no longer in trunk
	v1.52.1 // no longer in trunk
)

require github.com/stretchr/testify v1.8.0

require (
	github.com/davecgh/go-spew v1.1.1 // indirect
	github.com/pmezard/go-difflib v1.0.0 // indirect
	gopkg.in/yaml.v3 v3.0.1 // indirect
)
