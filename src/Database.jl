# using DBInterface
# using SQLite

# function initialize_database(sectors_filepath::AbstractString, temp_filepath::AbstractString=":memory:")
#     conn = DBInterface.connect(SQLite.DB, sectors_filepath)
#     SQLite.busy_timeout(conn, 20000)
#     DBInterface.execute(conn, "PRAGMA journal_mode=DELETE")
#     if DBInterface.execute(conn, "SELECT * from sqlite_master WHERE name=?", ("sectors",)) |> isempty
#         throw(ArgumentError("table sectors missing"))
#     elseif DBInterface.execute(conn, "SELECT * from sqlite_master WHERE name=?", ("schedule",)) |> isempty
#         throw(ArgumentError("table schedule missing"))
#     end

#     SQL(args...; kwargs...) = DBInterface.execute(conn, args...; kwargs...)

#     SQL("ATTACH DATABASE ? AS temporary", (temp_filepath,))

#     for db in ["", "temporary."]

#         if isempty(SQL("SELECT * FROM $(db)sqlite_master WHERE name=?", ("dense_energy_shifts",)))
#             DBInterface.execute(conn, """
#                 CREATE TABLE $(db)dense_energy_shifts (
#                     idx INTEGER NOT NULL,
#                     hopping REAL NOT NULL,
#                     interaction REAL NOT NULL,
#                     base_energy REAL NOT NULL,
#                     UNIQUE(idx, hopping, interaction),
#                     FOREIGN KEY(idx) REFERENCES sectors(idx)
#                 )
#             """)
#         end
#         if isempty(SQL("SELECT * FROM $(db)sqlite_master WHERE name=?", ("dense_eigen_results",)))
#             DBInterface.execute(conn, """
#                 CREATE TABLE $(db)dense_eigen_results (
#                     idx INTEGER NOT NULL,
#                     hopping REAL NOT NULL,
#                     interaction REAL NOT NULL,
#                     eigenindex INTEGER NOT NULL,
#                     eigenvalue REAL NOT NULL,
#                     UNIQUE(idx, hopping, interaction, eigenindex),
#                     FOREIGN KEY(idx) REFERENCES sectors(idx)
#                 )
#             """)
#         end
#         if isempty(SQL("SELECT * FROM $(db)sqlite_master WHERE name=?", ("dense_results",)))
#             DBInterface.execute(conn, """
#                 CREATE TABLE $(db)dense_results (
#                     idx INTEGER NOT NULL,
#                     hopping REAL NOT NULL,
#                     interaction REAL NOT NULL,
#                     temperature REAL NOT NULL,
#                     partition REAL NOT NULL,
#                     energy REAL NOT NULL,
#                     energy_squared REAL NOT NULL,
#                     correlation TEXT,
#                     UNIQUE(idx, hopping, interaction, temperature),
#                     FOREIGN KEY(idx) REFERENCES sectors(idx)
#                 )
#             """)
#         end

#         if isempty(SQL("SELECT * FROM $(db)sqlite_master WHERE name=?", ("sparse_results",)))
#             SQL("""
#                 CREATE TABLE IF NOT EXISTS $(db)sparse_results (
#                     idx INTEGER NOT NULL,
#                     hopping REAL NOT NULL,
#                     interaction REAL NOT NULL,
#                     temperature REAL NOT NULL,
#                     partition REAL NOT NULL,
#                     energy REAL NOT NULL,
#                     energy_squared REAL NOT NULL,
#                     correlation TEXT,
#                     samplecount INTEGER NOT NULL,
#                     UNIQUE(idx, hopping, interaction, temperature)
#                 )"""
#             )
#         end

#         if isempty(SQL("SELECT * FROM $(db)sqlite_master WHERE name=?", ("sparse_eigen_results",)))
#             DBInterface.execute(conn, """
#                 CREATE TABLE $(db)sparse_eigen_results (
#                     idx INTEGER NOT NULL,
#                     hopping REAL NOT NULL,
#                     interaction REAL NOT NULL,
#                     eigenindex INTEGER NOT NULL,
#                     eigenvalue REAL NOT NULL,
#                     run INTEGER NOT NULL,
#                     FOREIGN KEY(idx) REFERENCES sectors(idx)
#                 )
#             """)
#         end

#         if isempty(SQL("SELECT * FROM $(db)sqlite_master WHERE name=?", ("sparse_energy_shifts",)))
#             SQL("""
#                 CREATE TABLE IF NOT EXISTS sparse_energy_shifts (
#                     idx INTEGER NOT NULL,
#                     hopping REAL NOT NULL,
#                     interaction REAL NOT NULL,
#                     base_energy REAL NOT NULL,
#                     UNIQUE(idx, hopping, interaction),
#                     FOREIGN KEY(idx) REFERENCES sectors(idx)
#                 )
#             """)
#         end

#         #=
#         if isempty(SQL("SELECT * FROM $(db)sqlite_master WHERE name=?", ("sparse_results",)))
#             DBInterface.execute(conn, """
#                 CREATE TABLE $(db)sparse_results (
#                     idx INTEGER NOT NULL,
#                     hopping REAL NOT NULL,
#                     interaction REAL NOT NULL,
#                     seed INTEGER NOT NULL,
#                     samplecount INTEGER NOT NULL,
#                     filename TEXT NOT NULL,
#                     UNIQUE(idx, hopping, interaction, seed),
#                     FOREIGN KEY(idx) REFERENCES sectors(idx)
#                 )
#             """)
#         end
#         =#
#     end

#     return conn
# end

